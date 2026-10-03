///////////////////////////////////////////////////////////////////////////////////////////
// mml_export.h - Cross-platform symbol visibility macros for MML Packages
///////////////////////////////////////////////////////////////////////////////////////////
//
// Usage in public headers:
//   class MML_API MyClass { ... };
//   MML_API void myFunction();
//
// Per-package variants (for separate DLL builds):
//   class MML_FOURIER_API FourierTransform { ... };
//   class MML_OPTIMIZATION_API GeneticAlgorithm { ... };
//
// Build modes:
//   - Static library:  Define MML_STATIC (no import/export needed)
//   - Shared library:  MML_BUILDING_<PKG> set by CMake on library target
//   - Header-only:     Define MML_HEADER_ONLY (backward compatibility)
//
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_EXPORT_H
#define MML_EXPORT_H

// =============================================================================
// Header-only backward compatibility
// =============================================================================
#ifdef MML_HEADER_ONLY
    #define MML_API
    #define MML_FOURIER_API
    #define MML_OPTIMIZATION_API
    #define MML_PDE_API
    #define MML_STATISTICS_API
    #define MML_SYMBOLIC_API

// =============================================================================
// Static library build (no import/export decorations needed)
// =============================================================================
#elif defined(MML_STATIC)
    #define MML_API
    #define MML_FOURIER_API
    #define MML_OPTIMIZATION_API
    #define MML_PDE_API
    #define MML_STATISTICS_API
    #define MML_SYMBOLIC_API

// =============================================================================
// Shared library build (platform-specific visibility)
// =============================================================================
#else

    // -------------------------------------------------------------------------
    // MSVC: __declspec(dllexport/dllimport)
    // -------------------------------------------------------------------------
    #if defined(_MSC_VER)
        #ifdef MML_BUILDING_ALL
            #define MML_API __declspec(dllexport)
        #else
            #define MML_API __declspec(dllimport)
        #endif

        #ifdef MML_BUILDING_FOURIER
            #define MML_FOURIER_API __declspec(dllexport)
        #else
            #define MML_FOURIER_API __declspec(dllimport)
        #endif

        #ifdef MML_BUILDING_OPTIMIZATION
            #define MML_OPTIMIZATION_API __declspec(dllexport)
        #else
            #define MML_OPTIMIZATION_API __declspec(dllimport)
        #endif

        #ifdef MML_BUILDING_PDE
            #define MML_PDE_API __declspec(dllexport)
        #else
            #define MML_PDE_API __declspec(dllimport)
        #endif

        #ifdef MML_BUILDING_STATISTICS
            #define MML_STATISTICS_API __declspec(dllexport)
        #else
            #define MML_STATISTICS_API __declspec(dllimport)
        #endif

        #ifdef MML_BUILDING_SYMBOLIC
            #define MML_SYMBOLIC_API __declspec(dllexport)
        #else
            #define MML_SYMBOLIC_API __declspec(dllimport)
        #endif

    // -------------------------------------------------------------------------
    // GCC/Clang: __attribute__((visibility))
    // -------------------------------------------------------------------------
    #elif defined(__GNUC__) || defined(__clang__)
        #define MML_API __attribute__((visibility("default")))
        #define MML_FOURIER_API __attribute__((visibility("default")))
        #define MML_OPTIMIZATION_API __attribute__((visibility("default")))
        #define MML_PDE_API __attribute__((visibility("default")))
        #define MML_STATISTICS_API __attribute__((visibility("default")))
        #define MML_SYMBOLIC_API __attribute__((visibility("default")))

    // -------------------------------------------------------------------------
    // Unknown compiler: no decorations
    // -------------------------------------------------------------------------
    #else
        #define MML_API
        #define MML_FOURIER_API
        #define MML_OPTIMIZATION_API
        #define MML_PDE_API
        #define MML_STATISTICS_API
        #define MML_SYMBOLIC_API
    #endif

#endif // MML_STATIC / MML_HEADER_ONLY

#if defined(_MSC_VER)
    // Windows exports compiled template members via WINDOWS_EXPORT_ALL_SYMBOLS.
    #define MML_PDE_TEMPLATE_API
#else
    #define MML_PDE_TEMPLATE_API MML_PDE_API
#endif

#endif // MML_EXPORT_H