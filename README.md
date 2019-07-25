# AutoD

Calculate the value of function, at the same time, calculating the (arbitrary higher order) derivative.

# Header Only Project
This is a header only project.
To use this library. just add
```
#include "AutoD.h"
```
in front of your source file.

# Example
```
#include "AutoD.h"
int main() {
using namespace ad;        // The namespace of the library.
Real<2> x = Primary<2>(1); // Primary is a `making` function.
                           // Real<2> store the value and the first order, now it is (x, dx/dx, d^x/dx^2) = (1, 1, 0)

Real<2> y = 1;             // y is (1, 0, 0)
                           // y is not variable
PrimaryReset(y, 1);        // now y store the (0, 1, 0)
                           // the y is now variable

Real<2> xplusy = x + y;    // xplus y is (2, 2, 0)
                           // All primary variables represent one variable.
                           // It is by nature of light weight design of this library.

}
```

Currently, we supoort
```
+
-
*
/
sqrt
log
POW2
POW4
```
This project intents to show the potential usage auto-derivative calculation.


