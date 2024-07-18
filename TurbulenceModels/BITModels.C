// ADD THE MODELS TO THE RUNTIME TABLE //

#include "phaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// Typedefs
defineTurbulenceModelTypes
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    phaseModel
);

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (phaseModelPhaseCompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (phaseModelPhaseCompressibleTurbulenceModel, LES, Type)

#include "RAS/kEpsilonMa/kEpsilonMa.H"
makeRASModel(kEpsilonMa);

#include "RAS/kEpsilonMagolan/kEpsilonMagolan.H"
makeRASModel(kEpsilonMagolan);