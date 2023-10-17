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

Class
    Foam::wallLubricationModels::Hosokawa

Description
    Wall lubrication model of Hosokawa et al. (2002).

    References:
    \verbatim
        Hosokawa, S., Tomiyama, A., Misaki, S., & Hamada, T. (2002, January).
        Lateral migration of single bubbles due to the presence of wall.
        In ASME 2002 Joint US-European Fluids Engineering Division Conference
        (pp. 855-860).
        American Society of Mechanical Engineers.
    \endverbatim

SourceFiles
    Hosokawa.C

\*---------------------------------------------------------------------------*/

#ifndef Hosokawa_H
#define Hosokawa_H

#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                            Class Hosokawa Declaration
\*---------------------------------------------------------------------------*/

class Hosokawa
:
    public wallLubricationModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;


public:

    //- Runtime type information
    TypeName("Hosokawa");


    // Constructors

        //- Construct from components
        Hosokawa
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Hosokawa();


    // Member Functions

        //- Return phase-intensive wall lubrication force
        tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
