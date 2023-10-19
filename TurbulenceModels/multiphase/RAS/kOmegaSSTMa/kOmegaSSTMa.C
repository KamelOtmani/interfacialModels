/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include "kOmegaSSTMa.H"
#include "fvOptions.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTMa<BasicTurbulenceModel>::kOmegaSSTMa
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
    kOmegaSST<BasicTurbulenceModel>
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
bool kOmegaSSTMa<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
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
kOmegaSSTMa<BasicTurbulenceModel>::gasTurbulence() const
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
void kOmegaSSTMa<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2
)
{
    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    volScalarField yPlus
    (
        pow(this->betaStar_, 0.25)*this->y_*sqrt(this->k_)/this->nu()
    );

    this->nut_ =
        this->a1_*this->k_
       /max
        (
            this->a1_*this->omega_,
            this->b1_*this->F23()*sqrt(S2)
        )
      + sqr(1 - exp(-yPlus/16.0))
       *Cmub_*gasTurbulence.transport().d()*gasTurbulence.alpha()
       *(mag(this->U_ - gasTurbulence.U()));

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template <class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTMa<BasicTurbulenceModel>::kSource() const
{
    // return Ck()*dragModel_.K()*magSqr(pair_.Ur());

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

template <class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTMa<BasicTurbulenceModel>::omegaSource() const
{
    // return Ceps * k;
    // Ceps = 
    // return
    //     0.3
    //    *dragModel_.CdRe()
    //    *pair_.continuous().thermo().nu()
    //    /sqr(pair_.dispersed().d());

    
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const PhaseCompressibleTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());

    return
        alpha*rho*this->C3_*this->epsilon_*bubbleG()/this->k_
      + phaseTransferCoeff*gasTurbulence.epsilon()
      - fvm::Sp(phaseTransferCoeff, this->epsilon_);
}

template<class BasicTurbulenceModel>
void kOmegaSSTMa<BasicTurbulenceModel>::correct()
{
    kOmegaSST<BasicTurbulenceModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
