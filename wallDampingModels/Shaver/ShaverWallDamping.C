/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "ShaverWallDamping.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallDampingModels
{
    defineTypeNameAndDebug(Shaver, 0);
    addToRunTimeSelectionTable
    (
        wallDampingModel,
        Shaver,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::wallDampingModels::Shaver::limiter() const
{
    auto y = yWall();
    auto db = Cd_*pair_.dispersed().d();
    auto ydb = y/db;
    Foam::tmp<Foam::volScalarField> result;

    Foam::tmp<Foam::volScalarField>
    (
        new Foam::volScalarField
        (
            IOobject
            (
                "res",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("res",dimensionSet(0,0,0,0,0,0,0),0)
        )
    );

    forAll(ydb,faceI)
    {
        if ( ydb[faceI] < 0.5 ) {
            result[faceI] =  yWall() + scalar(0);
        }
        else if ( ydb > 1 ) {
            result[faceI] =  yWall() + scalar(1);
        }
        else {
            result[faceI] = (3 * pow(2 * ydb -1,2) - 2 * pow(2 * ydb -1,3));
        }
    }
    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDampingModels::Shaver::Shaver
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interpolated(dict, pair),
    Cd_("Cd", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDampingModels::Shaver::~Shaver()
{}


// ************************************************************************* //
