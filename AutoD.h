#pragma once

#include <stdint.h>
#include <vector>
#include <typeinfo>
#include <type_traits>
#include <cmath>

namespace ad
{

	template<int N>
	struct Real
	{
	public:
		//static_assert(N >= 0, "");
		Real() : Real(0) { }
		Real(double v) {
			fV[0] = v;
			for (int i = 1; i < Real<N>::fN; ++i) {
				fV[i] = 0;
			}
		}
		double V() { return fV[0]; }
		double D(int i) { return fV[i]; }

		//	private:
		static const int fN = N + 1;
		double fV[fN];
	};

	template<int N>
	Real<N> Primary(double v) {
		Real<N> res;
		if (Real<N>::fN >= 1) res.fV[0] = v;
		if (Real<N>::fN >= 2) res.fV[1] = 1;
		for (int i = 2; i < Real<N>::fN; ++i) {
			res.fV[i] = 0;
		}
		return res;
	}

	template<int N>
	void PrimanryReset(Real<N> &r, double v) {
		r = Primary<N>(v);
	}

	namespace dts {
		struct Zero
		{
			typedef Zero Diff;
			Zero diff() const { return Zero(); }
			double v() const { return 0; }

		};
		struct One
		{
			One() { }
			typedef Zero Diff;
			Zero diff() const { return Zero(); }
			double v() const { return 1; }

		};

		struct Const
		{
			double const fV;
			Const(double c) : fV(c) { }

			typedef Zero Diff;
			Diff diff() const { return Zero(); }
			double v() const { return fV; }

		};


		template<int n, int N>
		struct RealDiff {
			Real<N> const &v1_;
			RealDiff(Real<N> const &v) : v1_(v) { }

			typedef RealDiff<n + 1, N> Diff;

			auto diff() const {
				return RealDiff<n + 1, N>(v1_);
			}

			double v() const {
				return v1_.fV[n];
			}
		};

		template<class T1, class T2>
		struct Add;

		template<class T1, class T2>
		struct AddTraits {
			typedef Add<T1, T2> Type;
			static Type v(T1 const &v1, T2 const &v2) { return Type(v1, v2); }
		};
		template<class T1>
		struct AddTraits<T1, Zero> {
			typedef T1 Type;
			static Type v(T1 const &v1, Zero const &) { return v1; }
		};
		template<class T2>
		struct AddTraits<Zero, T2> {
			typedef T2 Type;
			static Type v(Zero const &, T2 const &v2) { return v2; }
		};

		template<class T1, class T2>
		struct Add {


			T1 const add1;
			T2 const add2;
			double const v_;
			typedef typename T1::Diff Diff1;
			typedef typename T2::Diff Diff2;
			typedef typename AddTraits<Diff1, Diff2>::Type Diff;

			auto diff() const {
				return AddTraits<Diff1, Diff2>::v(add1.diff(), add2.diff());
			}
			Add(T1 const &v1, T2 const &v2) : add1(v1), add2(v2), v_(add1.v() + add2.v()) {
			}
			double v() const {
				return v_;
			}
		};

		template<class T1, class T2>
		struct Sub;

		template<class T1, class T2>
		struct SubTraits {
			typedef Sub<T1, T2> Type;
			static Type v(T1 const &v1, T2 const &v2) { return Type(v1, v2); }
		};

		template<class T1>
		struct SubTraits<T1, Zero> {
			typedef T1 Type;
			static Type v(T1 const &v1, Zero const &) { return v1; }
		};

		template<class T1, class T2>
		struct Sub {

			T1 const sub1;
			T2 const sub2;
			double const v_;
			Sub(T1 const &v1, T2 const &v2) : sub1(v1), sub2(v2), v_(sub1.v() - sub2.v()) {
			}
			typedef typename T1::Diff Diff1;
			typedef typename T2::Diff Diff2;
			typedef typename SubTraits<Diff1, Diff2>::Type Diff;

			Diff diff() const {
				return SubTraits<Diff1, Diff2>::v(sub1.diff(), sub2.diff());
			}

			double v() const {
				return v_;
			}

		};

		template<class T1, class T2>
		struct Times;

		template<class T1, class T2>
		struct TimesTraits {
			typedef Times<T1, T2> Type;
			static Type v(T1 const &v1, T2 const &v2) { return Times<T1, T2>(v1, v2); }
		};

		template<class T1>
		struct TimesTraits<T1, Zero> {
			typedef Zero Type;
			static Type v(T1 const &, Zero const &) { return Zero(); }
		};
		template<class T1>
		struct TimesTraits<Zero, T1> {
			typedef Zero Type;
			static Type v(Zero const &, T1 const &) { return Zero(); }
		};

		template<class T1>
		struct TimesTraits<T1, One> {
			typedef T1 Type;
			static Type v(T1 const &v1, One const &) { return v1; }
		};

		template<class T2>
		struct TimesTraits<One, T2> {
			typedef T2 Type;
			static Type v(One const &, T2 const & v2) { return v2; }
		};

		template<class T1, class T2>
		struct Times {
			T1 const times1;
			T2 const times2;
			double const v_;

			typedef typename T1::Diff Diff1;
			typedef typename T2::Diff Diff2;
			typedef AddTraits<typename TimesTraits<Diff1, T2>::Type, typename TimesTraits<T1, Diff2>::Type > DiffTraits;
			typedef typename DiffTraits::Type Diff;

			auto diff() const {
				return DiffTraits::v(TimesTraits<Diff1, T2>::v(times1.diff(), times2),
					TimesTraits<T1, Diff2>::v(times1, times2.diff())
					);
			}

			Times(T1 const &v1, T2 const &v2) : times1(v1), times2(v2), v_(times1.v() * times2.v()) {
			}

			double v() const {
				return v_;
			}
		};


		template<class T1, int n>
		struct PowN;

		template<class T1, int n>
		struct PowNTraits {
			typedef PowN<T1, n> Type;
			static Type v(T1 const &v1) { return Type(v1); }
		};

		template<class T1>
		struct PowNTraits<T1, 1> {
			typedef T1 Type;
			static Type v(T1 const &v1) { return v1; }
		};

		template<class T1>
		struct PowNTraits<T1, 0> {
			typedef One Type;
			static Type v(T1 const &) { return One(); }
		};

		template<class T1, int n>
		struct PowN {
			T1 const pown;
			double const v_;
			PowN(T1 const &v1) : pown(v1), v_(cal_v()) { }

			typedef typename T1::Diff Diff1;
			typedef typename PowNTraits<T1, n - 1>::Type PowNN1;
			typedef TimesTraits<Const,
				typename TimesTraits<PowNN1, Diff1>::Type
			> DiffTraits;
			typedef typename DiffTraits::Type Diff;

			Diff diff() const {
				return DiffTraits::v(Const(n),
					TimesTraits<PowNN1, Diff1>::v(PowNTraits<T1, n - 1>::v(pown), pown.diff())
					);
			}

			double v() const {
				return v_;
			}

			double cal_v() const {
				if (n == 0) return 1;
				else if (n == 1) return pown.v();
				else if (n == 2) return pown.v()*pown.v();
				else if (n == 3) return pown.v()*pown.v()*pown.v();
				else if (n == -1) return 1 / pown.v();
				else if (n == -2) return 1 / (pown.v()*pown.v());
				else if (n == -3) return 1 / (pown.v()*pown.v()*pown.v());
				return pow(pown.v(), n);
			}

		};



		template<class T1, class T2>
		struct Div {
			T1 const div1;
			T2 const div2;
			double const v_;
			Div(T1 const &v1, T2 const &v2) : div1(v1), div2(v2), v_(div1.v()/div2.v()){
			}

			typedef typename T1::Diff Diff1;
			typedef typename T2::Diff Diff2;
			typedef typename TimesTraits<T1, Diff2>::Type MNType;

			typedef SubTraits<
				Div<Diff1, T2>,
				typename TimesTraits<MNType, PowN<T2, -2> >::Type
			> DiffTraits;
			typedef typename DiffTraits::Type Diff;

			Diff diff() const {
				return DiffTraits::v(
					Div<Diff1, T2>(div1.diff(), div2),
					Times<MNType, PowN<T2, -2> >(TimesTraits<T1, Diff2>::v(div1, div2.diff()), PowN<T2, -2>(div2))
				);
			}

			double v() const {
				return v_;
			}

		};

		template<class T1>
		struct Sqrt
		{

#if 1

			T1 const arg_;
			Sqrt(T1 const &v1) : arg_(v1) {
			}
			double v() const {
				return sqrt(arg_.v());
			}

#else
			T1 const arg_;
			double const v_;
			Sqrt(T1 const &v1) : arg_(v1), v_(sqrt(arg_.v())) {
			}
			double v() const {
				return v_;
			}
#endif


			typedef typename T1::Diff Diff1;
			typedef Div< Times<Const, Diff1>, Sqrt<T1> > Diff;

			Diff diff() const {
				return Diff(Times<Const, Diff1>(Const(0.5), arg_.diff()), *this);
			}


		};

		template<class T1>
		struct Log
		{
			T1 const arg_;
#if 0
			double const v_;
			Log(T1 const &v1) : arg_(v1), v_(log(arg_.v())) {
			}
			double v() const {
				return v_;
			}
#else
			Log(T1 const &v1) : arg_(v1) {
			}
			double v() const {
				return log(arg_.v());
			}
#endif

			typedef Div< typename T1::Diff, T1 > Diff;

			Diff diff() const {
				return Diff(arg_.diff(), arg_);
			}


		};

		template<int I, int E>
		struct Sm {
			template<int N, class Expr>
			static void v(Real<N> &res, Expr &m) {
				auto df = m.diff();
				res.fV[I] = df.v();
				Sm<I + 1, E>::v(res, df);
			}
		};

		template<int E>
		struct Sm<E, E> {
			template<int N, class Expr>
			static void v(Real<N> &, Expr  &) {
			}
		};

		template<int N, class T>
		Real<N> ToReal(T &v) {
			Real<N> res;
			res.fV[0] = v.v();
			Sm<1, Real<N>::fN>::v(res, v);
			return res;
		}

	}

	template<int N>
	Real<N> operator-(Real<N> const &r) {
		return Real<N>(0) - r;
	}

	template<int N>
	Real<N> operator+(Real<N> const &l, Real<N> const &r) {
		Real<N> res;
		for (int i = 0; i < Real<N>::fN; ++i) {
			res.fV[i] = l.fV[i] + r.fV[i];
		}
		return res;
	}

	template<int N>
	Real<N> operator+(double l, Real<N> const &r) {
		return Real<N>(l) + r;
	}

	template<int N>
	Real<N> operator+(Real<N> const &l, double r) {
		return l + Real<N>(r);
	}

	template<int N>
	Real<N> operator-(Real<N> const &l, Real<N> const &r) {
		Real<N> res;
		for (int i = 0; i < Real<N>::fN; ++i) {
			res.fV[i] = l.fV[i] - r.fV[i];
		}
		return res;
	}

	template<int N>
	Real<N> operator-(double l, Real<N> const &r) {
		return Real<N>(l) - r;
	}

	template<int N>
	Real<N> operator-(Real<N> const &l, double r) {
		return l - Real<N>(r);
	}

	template<int N>
	Real<N> operator*(Real<N> const &l, Real<N> const &r)
	{
		using namespace dts;
		Times<RealDiff<0, N>, RealDiff<0, N> > mul(l, r);
		return ToReal<N>(mul);
	}

	template<int N>
	Real<N> operator*(double l, Real<N> const &r) {
		return Real<N>(l) * r;
	}

	template<int N>
	Real<N> operator*(Real<N> const &l, double r) {
		return l * Real<N>(r);
	}

	template<int N>
	Real<N> operator/(Real<N> const &l, Real<N> const &r)
	{
		using namespace dts;
		Div<RealDiff<0, N>, RealDiff<0, N> > div(l, r);
		return ToReal<N>(div);
	}


	template<int N>
	Real<N> operator/(double l, Real<N> const &r) {
		return Real<N>(l) / r;
	}

	template<int N>
	Real<N> operator/(Real<N> const &l, double r) {
		return l / Real<N>(r);
	}

	template<int N>
	Real<N> sqrt(Real<N> const &l)
	{
		using namespace dts;
		Sqrt<RealDiff<0, N> > sqrt_(l);
		return ToReal<N>(sqrt_);
	}

	template<int N>
	Real<N> log(Real<N> const &l)
	{
		using namespace dts;
		Log<RealDiff<0, N> > log_(l);
		return ToReal<N>(log_);
	}

	template<int i>
	inline Real<i> POW2(Real<i> const &x) { return x*x; }
	template<int i>
	inline Real<i> POW4(Real<i> const &x) { return x*x*x*x; }

}
