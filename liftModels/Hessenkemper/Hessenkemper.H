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

Class
    Foam::liftModels::Hessenkemper

Description
    Lift model of Hessenkemper et al. (2021).

    Note that the dimensionless shear rate Sr is computed as described in
    Shi and Rzehak (2019).

    References:
    \verbatim
        Hessenkemper, H., Ziegenhein, T., Rzehak, R., Lucas, D.,
        & Tomiyama, A. (2021).
        Lift force coefficient of ellipsoidal single bubbles in water.
        International Journal of Multiphase Flow, 138, 103587.
    \endverbatim

    \verbatim
        Shi, P., & Rzehak, R. (2019).
        Lift forces on Solid Spherical Particles in Unbounded Flows.
        Chemical Engineering Science, 208, 115145.
    \endverbatim

SourceFiles
    Hessenkemper.C

\*---------------------------------------------------------------------------*/

#ifndef Hessenkemper_H
#define Hessenkemper_H

#include "liftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;
class aspectRatioModel;

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                        Class Hessenkemper Declaration
\*---------------------------------------------------------------------------*/

class Hessenkemper
:
    public liftModel
{
    // Private Data

        //- The aspect ratio model
        autoPtr<aspectRatioModel> aspectRatio_;

        //- Residual Reynolds number
        const dimensionedScalar residualRe_;


public:

    //- Runtime type information
    TypeName("Hessenkemper");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Hessenkemper
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Hessenkemper();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
