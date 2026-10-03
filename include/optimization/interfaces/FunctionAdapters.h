///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionAdapters.h                                                  ///
///  Description: Adapters bridging MML fixed-dimension functions into the dynamic    ///
///               (runtime-dimension) optimization objective interface.               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_OPTIMIZATION_FUNCTION_ADAPTERS_H
#define MML_OPTIMIZATION_FUNCTION_ADAPTERS_H

#include "MMLBase.h"
#include "base/Vector/VectorN.h"
#include "interfaces/IFunction.h"

#include <functional>

#include "Variables.h"   // OptimizationProblemError

namespace MML::Optimization {

	///////////////////////////////////////////////////////////////////////////
	///              FIXED-N  ->  DYNAMIC OBJECTIVE ADAPTER                  ///
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Bridge a fixed-dimension MML scalar function into a dynamic objective.
	 *
	 * The dynamic optimizers in this package take their objective as
	 * `std::function<Real(const Vector<Real>&)>`. To optimize an existing
	 * `IScalarFunction<N>` (compile-time dimension), wrap it with `AsDynamic`.
	 *
	 * `N` appears only here, at the call boundary; the optimizer core stays fully
	 * dynamic. The per-evaluation cost is a length-n copy, negligible relative to a
	 * typical objective evaluation. Because every algorithm in this package is
	 * derivative-free, nothing else drags `N` back into the framework.
	 *
	 * @tparam N Compile-time dimension of the wrapped function (its own property).
	 * @param f  The fixed-dimension function to adapt (must outlive the returned callable).
	 * @return A dynamic objective callable of runtime dimension N.
	 */
	template<int N>
	inline std::function<Real(const Vector<Real>&)> AsDynamic(const IScalarFunction<N>& f) {
		return [&f](const Vector<Real>& x) -> Real {
			if (static_cast<int>(x.size()) != N)
				throw OptimizationProblemError("AsDynamic: dimension mismatch");
			VectorN<Real, N> xn;
			for (int i = 0; i < N; ++i)
				xn[i] = x[i];
			return f(xn);
		};
	}

} // namespace MML::Optimization
#endif // MML_OPTIMIZATION_FUNCTION_ADAPTERS_H
