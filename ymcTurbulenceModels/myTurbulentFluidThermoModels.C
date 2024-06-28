/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "turbulentFluidThermoModels.H"

#include "CompressibleTurbulenceModel.H"
#include "CompressibleTurbulenceModel.H"
//#include "addToRuntimeSelectionTable.H"
//#include "makeTurbulenceModel.H"

//#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define createBaseTurbulenceModel(Alpha, Rho, baseModel, BaseModel, Transport) \
namespace Foam { \
typedef BaseModel<Transport> Transport##BaseModel; \
typedef RASModel<Transport> RAS##Transport##BaseModel; \
typedef LESModel<Transport> LES##Transport##BaseModel; \
}

 createBaseTurbulenceModel
 (
    geometricOneField,
    geometricOneField,
    compressibleTurbulenceModel,
    CompressibleTurbulenceModel,
    compressibleTransportModel
 );


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

  #include "iKOmegaSST.H"
  makeRASModel(iKOmegaSST);

// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //


// ************************************************************************* //
