/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "kEpsilonMagolan.H"
#include "fvOptions.H"
#include "twoPhaseSystem.H"
#include "dragModel.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonMagolan<BasicTurbulenceModel>::kEpsilonMagolan
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    gasTurbulencePtr_(nullptr),

    alphaInversion_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaInversion",
            this->coeffDict_,
            0.3
        )
    ),

    Cp_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cp",
            this->coeffDict_,
            0.25
        )
    ),

    C3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3",
            this->coeffDict_,
            this->C2_.value()
        )
    ),

    Cmub_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmub",
            this->coeffDict_,
            0.6
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonMagolan<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        alphaInversion_.readIfPresent(this->coeffDict());
        Cp_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        Cmub_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
const PhaseCompressibleTurbulenceModel
<
    typename BasicTurbulenceModel::transportModel
>&
kEpsilonMagolan<BasicTurbulenceModel>::gasTurbulence() const
{
    if (!gasTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& liquid = this->transport();
        const twoPhaseSystem& fluid =
            refCast<const twoPhaseSystem>(liquid.fluid());
        const transportModel& gas = fluid.otherPhase(liquid);

        gasTurbulencePtr_ =
           &U.db()
           .lookupObject<PhaseCompressibleTurbulenceModel<transportModel>>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    gas.name()
                )
            );
    }

    return *gasTurbulencePtr_;
}


template<class BasicTurbulenceModel>
void kEpsilonMagolan<BasicTurbulenceModel>::correctNut()
{
    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    //TODO: Probably store the ref to the gas phase so we dont need to keep copying this
    const transportModel& liquid = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(liquid.fluid());
    const transportModel& gas = fluid.otherPhase(liquid);

    // volScalarField nuBIT = Cmub_*gasTurbulence.transport().d()*gasTurbulence.alpha()
    //    *(mag(this->U_ - gasTurbulence.U()));
    // Cmu_B of Magolan 2018
    volScalarField CmuMagolan = 1.5 - 0.5 * exp(-10 * gas);
    this->nut_ =
        this->Cmu_*sqr(this->k_)/this->epsilon_ * CmuMagolan;

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonMagolan<BasicTurbulenceModel>::bubbleG() const
{
    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    const transportModel& liquid = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(liquid.fluid());
    const transportModel& gas = fluid.otherPhase(liquid);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    volScalarField magUr(mag(this->U_ - gasTurbulence.U()));
    volScalarField Re = (magUr*gas.d())/liquid.nu();
    volScalarField Cd = drag.CdRe()/(magUr*gas.d())*liquid.nu();
    // Lahey
    // volScalarField kBI = (4.0/(3.0*Cd))*(Cp_*(1+pow(Cd,4.0/3.0)));
    // Ma 2017
    // volScalarField kBI = min(0.18 * Re,scalar(1)) ;
    // Magolan 2019
    // volScalarField kBI = scalar(0.34);

    tmp<volScalarField> bubbleG
    (
        (3.0/4)*(
            pow3(magUr)
          * Cd
          * /*kBI = */0.34
        )
       *gas
       /gas.d()
    );
    return bubbleG;
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
kEpsilonMagolan<BasicTurbulenceModel>::phaseTransferCoeff() const
{
    const volVectorField& U = this->U_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const turbulenceModel& gasTurbulence = this->gasTurbulence();

    return
    (
        max(alphaInversion_ - alpha, scalar(0))
       *rho
       *min(gasTurbulence.epsilon()/gasTurbulence.k(), 1.0/U.time().deltaT())
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonMagolan<BasicTurbulenceModel>::kSource() const
{
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());

    return
        alpha*rho*bubbleG()
      + phaseTransferCoeff*gasTurbulence.k()
      - fvm::Sp(phaseTransferCoeff, this->k_);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonMagolan<BasicTurbulenceModel>::epsilonSource() const
{
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    const transportModel& liquid = this->transport();
    const twoPhaseSystem& fluid = refCast<const twoPhaseSystem>(liquid.fluid());
    const transportModel& gas = fluid.otherPhase(liquid);

    // const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    volScalarField magUr(mag(this->U_ - gasTurbulence.U()));

    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());
    // Lahey
    // volScalarField Tau = this->k_/this->epsilon_;
    // Ma
    // volScalarField Tau = gas.d()/magUr;
    // Magolan
    volScalarField Tau = gas.d()*pow(gas,1.0/3.0)/magUr;

    return
        alpha*rho*this->C3_*bubbleG()/Tau
      + phaseTransferCoeff*gasTurbulence.epsilon()
      - fvm::Sp(phaseTransferCoeff, this->epsilon_);
}


template<class BasicTurbulenceModel>
void kEpsilonMagolan<BasicTurbulenceModel>::correct()
{
    kEpsilon<BasicTurbulenceModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
