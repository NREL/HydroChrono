/*********************************************************************
 * @file  hydro_forces.h
 * @brief Backward-compatibility header - includes hydro_system.h
 *
 * This header is retained for backward compatibility. It simply includes
 * the new hydro_system.h header which defines HydroSystem (the user-facing
 * hydrodynamics façade) and the HydroForces alias.
 *
 * @deprecated Include <hydroc/hydro_system.h> directly in new code.
 *********************************************************************/

#ifndef HYDRO_FORCES_H
#define HYDRO_FORCES_H

// Forward to the new header
#include <hydroc/hydro_system.h>

#endif  // HYDRO_FORCES_H

