// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief As propriedades do problema em que H2 é injetado sob uma camada de baixa permeabilidade a uma profundidade de 2700m.
 */

#ifndef DUMUX_INJECTION_PROPERTIES_HH
#define DUMUX_INJECTION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/mpfa/omethod/staticinteractionvolume.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Criar novas tags de tipo
namespace TTag {
struct Injection { using InheritsFrom = std::tuple<TwoPTwoC>; };
struct InjectionBox { using InheritsFrom = std::tuple<Injection, BoxModel>; };
struct InjectionCCTpfa { using InheritsFrom = std::tuple<Injection, CCTpfaModel>; };
struct InjectionCCMpfa { using InheritsFrom = std::tuple<Injection, CCMpfaModel>; };
} // fim do namespace TTag

// Definir o tipo de malha
template<class TypeTag>
struct Grid<TypeTag, TTag::Injection> { using type = Dune::YaspGrid<2>; };

// Definir a propriedade do problema
template<class TypeTag>
struct Problem<TypeTag, TTag::Injection> { using type = InjectionProblem<TypeTag>; };

// Definir a configuração do fluido
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Injection>
{
    using type = FluidSystems::H2OH2<GetPropType<TypeTag, Properties::Scalar>,
                                     FluidSystems::H2OH2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

// Definir os parâmetros espaciais
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Injection>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InjectionSpatialParams<GridGeometry, Scalar>;
};

// Definir se frações molares (true) ou mássicas (false) são usadas
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Injection> { static constexpr bool value = true; };

// Ativar ou não o cache (soluções de referência criadas sem cache)
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Injection> { static constexpr bool value = ENABLECACHING; };

// usar o volume de interação estático ao redor de vértices internos no teste mpfa
template<class TypeTag>
struct PrimaryInteractionVolume<TypeTag, TTag::InjectionCCMpfa>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;

    // malha bidimensional estruturada
    static constexpr int numIvScvs = 4;
    static constexpr int numIvScvfs = 4;

    // usar as características padrão
    using Traits = CCMpfaODefaultStaticInteractionVolumeTraits< NodalIndexSet, Scalar, numIvScvs, numIvScvfs >;
public:
    using type = CCMpfaOStaticInteractionVolume< Traits >;
};

} // fim do namespace Dumux::Properties

#endif
