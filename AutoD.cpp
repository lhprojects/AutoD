#include "AutoD.h"
#include <chrono>
#define TPrintf printf
const double PI = 3.1415926535897932384626433;

static int n_failed = 0;
#define TEST_SAME(x, y) do { double x_ = (x); double y_ = (y);\
	if(isnan(x_) || fabs(x_ - y_) > abs(y_)*1E-3)  { n_failed++; TPrintf("FAILED: %s: %3d:"  #x " (%f) == " #y " (%f)\n", __func__, __LINE__, x_, y_); }   } while(0)
#define TEST_TRUE(x) do { double x_ = (x); if(!x_)  { n_failed++; TPrintf("FAILED: %s: %3d:"  #x "\n", __func__, __LINE__);}   } while(0)

int main() {

	using namespace ad;
	using namespace  dts;

	{
		Real<2> x = Primary<2>(1);
		Real<2> y = x + x;
		TEST_SAME(y.V(), 2);
		TEST_SAME(y.D(1), 2);
	}
	{
		Real<2> x = Primary<2>(1);
		Real<2> y = x - x;
		TEST_SAME(y.V(), 0);
		TEST_SAME(y.D(1), 0);
	}

	{
		using namespace ad;
		Real<2> x = Primary<2>(3);
		typedef Add<RealDiff<0, 2>, RealDiff<0, 2> > AddT;
		typedef Times<RealDiff<0, 2>, RealDiff<0, 2> > TimesT;
		typedef Div<RealDiff<0, 2>, RealDiff<0, 2> > DivT;
		typedef PowN<RealDiff<0, 2>, 2 > PowNT;
		AddT add(x, x);
		TimesT times(x, x);
		DivT div(x, x);
		PowNT pownt(x);
		//DivT::Diff d = div.diff();
		DivT::Diff::Diff d = div.diff().diff();
		printf("%s\n\n", typeid(d).name());

		printf("%s\n\n", typeid(pownt).name());
		PowNT::Diff pownt_d = pownt.diff();
		printf("%s\n\n", typeid(pownt_d).name());
		PowNT::Diff::Diff pownt_dd = pownt.diff().diff();
		printf("%s\n\n", typeid(pownt_dd).name());

		TEST_SAME(add.v(), 6);
		TEST_SAME(add.diff().v(), 2);
		TEST_SAME(times.v(), 9);
		TEST_SAME(times.diff().v(), 6);

	}

	{
		Real<2> x = Primary<2>(3);
		Real<2> y = x * x;
		Real<2> z = x * x * x;


		TEST_SAME(y.V(), 9);
		TEST_SAME(y.D(1), 6);
		TEST_SAME(y.D(2), 2);

		TEST_SAME(z.V(), 27);
		TEST_SAME(z.D(1), 27);
		TEST_SAME(z.D(2), 18);

	}

	{
		Real<2> x = Primary<2>(3);
		Real<2> y = x / x;


		TEST_SAME(y.V(), 1);
		TEST_SAME(y.D(1), 0);
		TEST_SAME(y.D(2), 0);

	}

	{
		Real<2> x = Primary<2>(3);
		Real<2> y = sqrt(x);


		TEST_SAME(y.V(), sqrt(3));
		TEST_SAME(y.D(1), 0.5 / sqrt(3.));
		TEST_SAME(y.D(2), -0.5*0.5 / sqrt(3) / sqrt(3) / sqrt(3));

	}

	{ // log		
		Real<2> x = Primary<2>(3);
		Real<2> y = log(x);

		TEST_SAME(y.V(), log(3.0));
		TEST_SAME(y.D(1), 1 / 3.);
		TEST_SAME(y.D(2), -1 / 9.);

		Log< RealDiff<0, 2> > log_(x);
		printf("%s\n\n", typeid(y).name());
		printf("%s\n\n", typeid(log_.diff()).name());
		printf("%s\n\n", typeid(log_.diff().diff()).name());

	}

	{
		Real<3> x = Primary<3>(3);
		Real<3> y = x * x;
		Real<3> z = y * y;


		TEST_SAME(z.V(), 3 * 3 * 3 * 3);
		TEST_SAME(z.D(1), 4 * 3 * 3 * 3);
		TEST_SAME(z.D(2), 4 * 3 * 3 * 3);
		TEST_SAME(z.D(3), 4 * 3 * 2 * 3);

	}

	{
		Real<3> x = Primary<3>(rand());
		Real<3> y = x * x;
		Real<3> z = y * y;
		printf("%p\n", z.fV);


	}

	{
		auto t0 = std::chrono::high_resolution_clock::now();
		Real<4> res;
		for (int i = 0; i < 100000; ++i) {
			Real<4> x = Primary<4>(3);
			Real<4> x2 = x * x;
			Real<4> x2x = x * x + x;
			Real<4> x2xx = x * x - x;
			Real<4> x2xxsqrt = sqrt(x2xx);
			res = res + log(x2xxsqrt);
		}
		auto t1 = std::chrono::high_resolution_clock::now();
		printf("time used: %d us\n", (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
		printf("%f %f %f %f\n", res.V(), res.D(1), res.D(2), res.D(3));
	}


	printf("%d test(s) failed\n", n_failed);
	return 0;
}
