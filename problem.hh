// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Problema onde H2 é injetado sob uma camada de baixa permeabilidade a uma profundidade de 2700m.
 */

#ifndef DUMUX_INJECTION_PROBLEM_HH
#define DUMUX_INJECTION_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/h2oh2.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Problema onde H2 é injetado sob uma camada de baixa permeabilidade a uma profundidade de 2700m.
 *
 */
template <class TypeTag>
class InjectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    // Define aliases para tipos usados
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    // Índices das variáveis primárias
    enum
    {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx
    };

    // Presença de fases
    enum { wPhaseOnly = Indices::firstPhaseOnly };

    // Índices das equações
    enum
    {
        contiH2OEqIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        contiH2EqIdx = Indices::conti0EqIdx + FluidSystem::H2Idx,
    };

    // Índices das fases
    enum
    {
        gasPhaseIdx = FluidSystem::H2Idx,
        H2OIdx = FluidSystem::H2OIdx,
        H2Idx = FluidSystem::H2Idx
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    // Propriedade que define se frações molares ou mássicas são usadas
    static constexpr bool useMoles = ModelTraits::useMoles();

public:
    InjectionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // [MODIFICAÇÃO 1/3] - Adicionado parâmetro de salinidade
        salinity_ = getParam<Scalar>("Problem.Salinity", 0.0); // Valor padrão 0 (água doce)

        nTemperature_       = getParam<int>("Problem.NTemperature");
        nPressure_          = getParam<int>("Problem.NPressure");
        pressureLow_        = getParam<Scalar>("Problem.PressureLow");
        pressureHigh_       = getParam<Scalar>("Problem.PressureHigh");
        temperatureLow_     = getParam<Scalar>("Problem.TemperatureLow");
        temperatureHigh_    = getParam<Scalar>("Problem.TemperatureHigh");
        depthBOR_           = getParam<Scalar>("Problem.DepthBOR");
        name_               = getParam<std::string>("Problem.Name");

        // inicializa as tabelas do sistema de fluidos
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        // imprime no console se estão sendo usadas frações molares ou mássicas
        if(useMoles)
            std::cout<<"problem uses mole-fractions"<<std::endl;
        else
            std::cout<<"problem uses mass-fractions"<<std::endl;
    }

    /*!
     * \brief Retorna o nome do problema
     *
     * Usado como prefixo para os arquivos gerados pela simulação.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Calcula a densidade da água salgada considerando efeitos de salinidade
     * \param temperature Temperatura [K]
     * \param pressure Pressão [Pa]
     * \note Usa o valor de salinidade armazenado no problema
     * \note Implementa uma correção linear de densidade: ρ_salgada = ρ_pura * (1 + 0.7*salinidade)
     *       Válido para salinidade < 0.25 (25% em fração mássica)
     */
    Scalar salineWaterDensity(Scalar temperature, Scalar pressure) const
    {
        const Scalar rhoPure = FluidSystem::H2O::liquidDensity(temperature, pressure);
        //Correção linear para salinidade
        return rhoPure * (1.0 + 0.7 * salinity_);
    }

    /*!
     * \brief Especifica quais condições de contorno devem ser
     *        usadas para cada equação em um segmento de contorno
     *
     * \param globalPos A posição global
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (globalPos[0] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Avalia as condições de contorno para um segmento de fronteira Dirichlet
     *
     * \param globalPos A posição global
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Avalia as condições de contorno para um segmento de fronteira Neumann
     *        em função da solução atual.
     *
     * \param element O elemento finito
     * \param fvGeometry A geometria de volume finito do elemento
     * \param elemVolVars Todas as variáveis de volume do elemento
     * \param elemFluxVarsCache Cache das variáveis de fluxo para todas as faces no stencil
     * \param scvf A face de subvolume de controle
     *
     * Esse método é usado quando a condição de Neumann depende da
     * solução e requer quantidades específicas do método totalmente implícito.
     * Os valores armazenam o fluxo de massa de cada fase normal à fronteira.
     * Valores negativos indicam entrada (inflow).
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& globalPos = scvf.ipGlobal();

        Scalar injectedPhaseMass = 1e-3;
        Scalar moleFracW = elemVolVars[scvf.insideScvIdx()].moleFraction(gasPhaseIdx, H2OIdx);
        if (globalPos[1] < 14 - eps_ && globalPos[1] > 6.5 - eps_)
        {
            values[contiH2EqIdx] = -(1-moleFracW)*injectedPhaseMass/FluidSystem::molarMass(H2Idx);
            values[contiH2OEqIdx] = -moleFracW*injectedPhaseMass/FluidSystem::molarMass(H2OIdx);
        }
        return values;
    }

    /*!
     * \brief Avalia os valores iniciais para um volume de controle.
     *
     * \param globalPos A posição global
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }

private:
    /*!
     * \brief Avalia os valores iniciais para um volume de controle.
     *
     * Método interno para condição inicial
     *
     * \param globalPos A posição global
     */
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);

        //Substituído por chamada à salineWaterDensity()
        Scalar densityW = salineWaterDensity(
            this->spatialParams().temperatureAtPos(globalPos),
            1e5
        );

        Scalar pl = 1e5 - densityW*this->spatialParams().gravity(globalPos)[1]*(depthBOR_ - globalPos[1]);
        Scalar moleFracLiquidH2 = pl*0.95/BinaryCoeff::H2O_H2::henry(this->spatialParams().temperatureAtPos(globalPos));
        Scalar moleFracLiquidH2O = 1.0 - moleFracLiquidH2;

        Scalar meanM =
            FluidSystem::molarMass(H2OIdx)*moleFracLiquidH2O +
            FluidSystem::molarMass(H2Idx)*moleFracLiquidH2;
        if(useMoles)
        {
            // formulação em fração molar
            priVars[switchIdx] = moleFracLiquidH2;
        }
        else
        {
            // formulação em fração mássica
            Scalar massFracLiquidH2 = moleFracLiquidH2*FluidSystem::molarMass(H2Idx)/meanM;
            priVars[switchIdx] = massFracLiquidH2;
        }
        priVars[pressureIdx] = pl;
        return priVars;
    }

    Scalar depthBOR_;
    static constexpr Scalar eps_ = 1e-6;

    int nTemperature_;
    int nPressure_;
    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    Scalar salinity_; //Fração mássica de sal (0 = água doce)

    std::string name_;
};

} // fim do namespace Dumux

#endif
