// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Os parâmetros espaciais do problema onde ar é injetado sob uma camada de baixa permeabilidade a uma profundidade de 2700m.
 */

#ifndef DUMUX_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Definição dos parâmetros espaciais para o problema de injeção
 *        que utiliza o modelo isotérmico bifásico bicomponente totalmente implícito.
 */
template<class GridGeometry, class Scalar>
class InjectionSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                       InjectionSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ThisType = InjectionSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Exporta o tipo usado para a permeabilidade
    using PermeabilityType = Scalar;

    InjectionSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , finePcKrSwCurve_("SpatialParams.FineMaterial")
    , coarsePcKrSwCurve_("SpatialParams.CoarseMaterial")
    {
        layerBottom_ = 22.5;

        // permeabilidades intrínsecas
        fineK_ = 1e-13;
        coarseK_ = 1e-12;

        // porosidades
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // temperatura
        temperature_ = getParam<Scalar>("SpatialParams.InitialTemperature");
    }

    /*!
     * \brief Retorna o tensor de permeabilidade intrínseca \f$[m^2]\f$
     *
     * \param globalPos A posição global
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Retorna a porosidade \f$[-]\f$
     *
     * \param globalPos A posição global
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return finePorosity_;
        return coarsePorosity_;
    }

    /*!
     * \brief Retorna o objeto de parâmetros para a relação entre pressão capilar /
     *        saturação da matriz do material
     *
     * \param globalPos A posição global
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return makeFluidMatrixInteraction(finePcKrSwCurve_);
        return makeFluidMatrixInteraction(coarsePcKrSwCurve_);
    }

    /*!
     * \brief Função para definir qual fase é considerada como a fase molhante.
     *
     * \param globalPos A posição do centro do elemento
     * \return O índice da fase molhante
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

    /*!
     * \brief Retorna a temperatura no domínio para a posição fornecida
     * \param globalPos A posição em coordenadas globais onde a temperatura deve ser especificada
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return temperature_;
    }

private:
    // Define se a região é de material mais fino com base na altura
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > layerBottom_; }

    Scalar fineK_;
    Scalar coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    Scalar temperature_;

    const PcKrSwCurve finePcKrSwCurve_;
    const PcKrSwCurve coarsePcKrSwCurve_;
};

} // fim do namespace Dumux

#endif
