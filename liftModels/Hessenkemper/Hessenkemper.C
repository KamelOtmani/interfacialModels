/*---------------------------------------------------------------------------*\
  ==  == ====== ====   ====    |
                    \\     ||  | HZDR Multiphase Addon for OpenFOAM
  ======   //   ||  || ===//   | Website: https://doi.org/10.14278/rodare.767
  ||  ||  //    ||  // || \\   | Copyright (C) 2019-2021 OpenFOAM Foundation
  ==  == ====== ====   ==  ==  | and Helmholtz-Zentrum Dresden-Rossendorf
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Hessenkemper.H"
#include "phasePair.H"
#include "aspectRatioModel.H"
#include "fvcCurl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(Hessenkemper, 0);
    addToRunTimeSelectionTable(liftModel, Hessenkemper, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::Hessenkemper::Hessenkemper
(
    const dictionary& dict,
    const phasePair& pair
)
:
    liftModel(dict, pair),
    aspectRatio_(aspectRatioModel::New(dict.subDict("aspectRatio"), pair)),
    residualRe_
    (
        dimensionedScalar::lookupOrDefault("residualRe", dict, dimless, 1e-1)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::Hessenkemper::~Hessenkemper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::Hessenkemper::Cl() const
{
    // const volScalarField EoH
    // (
    //     pair_.Eo(pair_.dispersed().d()/cbrt(aspectRatio_->E()))
    // );

    volScalarField EoH(pair_.EoH2());

    volScalarField Re(max(pair_.Re(), residualRe_));

    volScalarField G
    (
        0.0275*log(1 + exp(4*(EoH - 5.6))) - 0.14*(EoH - 5.2) - 0.44
    );

    volScalarField fEoH
    (
        log(1 + exp(-12*G))/12
    );

    volScalarField Sr
    (
        sqr(pair_.dispersed().d())
       /(
            Re
           *pair_.continuous().thermo().nu()
        )
       *mag(fvc::curl(pair_.continuous().U()))
    );

    volScalarField ClLowSqr
    (
        sqr(6*2.255)
       *sqr(Sr)
       /(
            pow4(constant::mathematical::pi)
           *Re
           *pow3(Sr + 0.2*Re)
        )
    );

    volScalarField ClHighSqr
    (
        sqr(0.5*(Re + 16)/(Re + 29))
    );

    volScalarField fSrRe
    (
        sqrt(ClLowSqr + ClHighSqr)
    );

    return fSrRe - fEoH;
}


// ************************************************************************* //
