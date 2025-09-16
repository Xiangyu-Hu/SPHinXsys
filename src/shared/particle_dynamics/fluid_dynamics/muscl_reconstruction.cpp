
 #include "muscl_reconstruction.hpp"

 namespace SPH {
 namespace fluid_dynamics {
 
 static const char* to_cstr(SlopeLimiter L)
 {
     switch (L) {
         case SlopeLimiter::None:    return "None";
         case SlopeLimiter::Minmod:  return "Minmod";
         case SlopeLimiter::MC:      return "MC";
         case SlopeLimiter::VanLeer: return "VanLeer";
         default:                    return "Unknown";
     }
 }
 
 // (Optional) expose a tiny helper for logging
 const char* limiter_name(SlopeLimiter L) { return to_cstr(L); }
 
 } // namespace fluid_dynamics
 } // namespace SPH
 