// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Teste para o problema onde H2 é injetado sob uma camada de baixa permeabilidade a uma profundidade de 2700m.
 */
#include <config.h>

#include <iostream>

#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/loadsolution.hh>

#include <dumux/experimental/assembly/multistagefvassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define o type tag para este problema
    using TypeTag = Properties::TTag::TYPETAG;

    // inicializa MPI e/ou backend de multithread, se necessário
    Dumux::initialize(argc, argv);

    // processa os argumentos de linha de comando e o arquivo de entrada
    Parameters::init(argc, argv);

    // tenta criar uma malha (a partir de um arquivo fornecido ou das configurações do input)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // computa usando a visão da malha (grid) folha
    const auto& leafGridView = gridManager.grid().leafGridView();

    // cria a geometria de volume finito da malha
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // define o problema (condições iniciais e de contorno)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // obtém os parâmetros do laço de tempo
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // verifica se a simulação está sendo reiniciada após interrupção
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // vetor solução
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    if (restartTime > 0)
    {
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        const auto fileName = getParam<std::string>("Restart.File");
        const auto pvName = createPVNameFunction<IOFields, PrimaryVariables, ModelTraits, FluidSystem>();
        loadSolution(x, fileName, pvName, *gridGeometry);
    }
    else
        problem->applyInitialSolution(x);
    auto xOld = x;

    // variáveis da malha
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // inicializa o módulo de saída VTK
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // adiciona campos específicos do modelo
    vtkWriter.write(restartTime);

    // instancia o laço de tempo
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    auto timeSteppingMethod = std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();

    // montador com laço de tempo para problema instacionário
    using Assembler = Experimental::MultiStageFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeSteppingMethod, xOld);

    // solver linear
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // solver não-linear
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    // integrador de tempo multietapas
    using TimeStepper = Experimental::MultiStageTimeStepper<NewtonSolver>;
    TimeStepper timeStepper(nonLinearSolver, timeSteppingMethod);

    // laço de tempo principal
    timeLoop->start(); do
    {
        // integração no tempo
        timeStepper.step(x, timeLoop->time(), timeLoop->timeStepSize());

        // atualiza a solução anterior
        xOld = x;
        gridVariables->advanceTimeStep();

        // avança o laço de tempo
        timeLoop->advanceTimeStep();

        // exibe estatísticas do passo atual
        timeLoop->reportTimeStep();

        // ajusta novo dt sugerido pelo solver de Newton
        timeLoop->setTimeStepSize(nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

        // escreve saída VTK
        vtkWriter.write(timeLoop->time());

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    if (leafGridView.comm().rank() == 0)
        Parameters::print();

    return 0;
} // fim da função main
