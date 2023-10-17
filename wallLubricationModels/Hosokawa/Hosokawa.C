/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "Hosokawa.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLubricationModels
{
    defineTypeNameAndDebug(Hosokawa, 0);
    addToRunTimeSelectionTable
    (
        wallLubricationModel,
        Hosokawa,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLubricationModels::Hosokawa::Hosokawa
(
    const dictionary& dict,
    const phasePair& pair
)
:
    wallLubricationModel(dict, pair),
    residualRe_
    (
        dimensionedScalar::lookupOrDefault("residualRe", dict, dimless, 1e-1)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLubricationModels::Hosokawa::~Hosokawa()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::wallLubricationModels::Hosokawa::Fi() const
{
    const volVectorField Ur(pair_.Ur());

    const volVectorField& n(nWall());
    const volScalarField& y(yWall());
    const volScalarField Eo(pair_.Eo());

    volScalarField Re(max(pair_.Re(), residualRe_));

    return zeroGradWalls
    (
        max
        (
            dimensionedScalar(inv(dimLength), Zero),
            max
            (
                0.01085*Eo,
                3.5/pow(Re, 1.9)
            )
           *pair_.dispersed().d()/sqr(y)
        )
       *pair_.continuous().rho()
       *magSqr(Ur - (Ur & n)*n)
       *n
    );
}


// ************************************************************************* //
