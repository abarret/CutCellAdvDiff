// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "LS/LSCartGridFunction.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/app_namespaces.h"
#include <ibamr/config.h>

#include <SAMRAI_config.h>

namespace LS
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

LSCartGridFunction::LSCartGridFunction(const string& object_name) : CartGridFunction(object_name)
{
    // intentionally blank
    return;
} // LSCartGridFunction

} // namespace LS
