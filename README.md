# lusolve

C++ port of ruby's `LUSolve` module, solves `A・x = b` for `x`, using LU decomposition.

# dependency

C++11 compiler.

# usage

```c++
#include "lusolve.hpp"
#include <algorithm>

int main()
{
    int const kSize = 3;
    
    std::vector<std::vector<double>> a;
    a.resize(kSize);
    std::fill_n(a.begin(), kSize, std::vector<double>(kSize));
    a[0][0] = 1;
    a[0][1] = 2;
    // setup matrix "a".
    // ...

    std::vector<double> b;
    b.resize(kSize);
    // setup vector "b".
    // ...

    std::vector<double> x;
    try {
        x = lusolve(a, b);
    } catch (std::runtime_error const&) {
        // in case singular matrix
    }

    return 0;
}
```

# license

Same as ruby's license.

* https://github.com/ruby/ruby/blob/2cd6800fd8437b1f862f3f5c44db877159271d17/COPYING