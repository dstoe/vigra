/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2014 by                                 */
/*               Ullrich Koethe,                                        */
/*               Esteban Pardo                                          */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#define VIGRA_CHECK_BOUNDS

#include <algorithm>
#include <vigra/unittest.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/simplex.hxx>


using namespace vigra;

template <class Point>
struct UnitMatrix {
    UnitMatrix()
    : index_(0)
    {}

    Point operator()()
    {
        Point vec(0);
        if (index_ < vec.size()) {
            vec[index_] = 1;
            index_++;
        }
        return vec;
    }

    int index_;
};

struct SimplexTest
{
    typedef vigra::TinyVector<double, 1> Vector1;
    typedef vigra::TinyVector<double, 2> Vector2;
    typedef vigra::TinyVector<double, 3> Vector3;
    typedef vigra::Simplex<Vector1> Simplex1;
    typedef vigra::Simplex<Vector2> Simplex2;
    typedef vigra::Simplex<Vector3> Simplex3;
    typedef vigra::Facet<Vector2> Facet2;
    typedef vigra::Facet<Vector3> Facet3;

    template <int N>
    vigra::Simplex<vigra::TinyVector<double, N> > makeSimplex()
    {
        typedef typename vigra::TinyVector<double, N> Vector;
        typename vigra::ArrayVector<Vector> vertices(N + 1, Vector());
        std::generate(vertices.begin(), vertices.end(), UnitMatrix<Vector>());

        vigra::Simplex<Vector> ret(vertices.begin(), vertices.end());
        shouldEqual(ret.size(), N + 1);
        return ret;
    }

    void testBasics()
    {
        // 1D
        Simplex1 simplex1(makeSimplex<1>());

        shouldEqual(simplex1.size(), 2);
        shouldEqual(simplex1.dimension, 1);
        shouldEqual(simplex1.nVolume(), 1.0);
        shouldEqual(simplex1.contains(Vector1(0.5)), true);     // inner
        shouldEqual(simplex1.contains(Vector1(1.0)), true);     // vertex
        shouldEqual(simplex1.contains(Vector1(1.5)), false);    // outer

        // 2D
        Simplex2 simplex2(makeSimplex<2>());

        shouldEqual(simplex2.size(), 3);
        shouldEqual(simplex2.dimension, 2);
        shouldEqual(simplex2.nVolume(), 0.5);
        shouldEqual(simplex2.contains(Vector2(0.2, 0.2)), true);    // inner
        shouldEqual(simplex2.contains(Vector2(1.0, 0.0)), true);    // vertex
        shouldEqual(simplex2.contains(Vector2(0.5, 0.5)), true);    // edge
        shouldEqual(simplex2.contains(Vector2(1.0, 1.0)), false);   // outer

        Facet2 facet2(simplex2, 0);
        shouldEqual(facet2.size(), 2);
        shouldEqual(facet2.nSurface(), 1.0);

        // 3D
        Simplex3 simplex3(makeSimplex<3>());

        shouldEqual(simplex3.size(), 4);
        shouldEqual(simplex3.dimension, 3);
        shouldEqual(simplex3.nVolume(), 1.0 / 6.0);
        shouldEqual(simplex3.contains(Vector3(0.1, 0.1, 0.1)), true);   // inner
        shouldEqual(simplex3.contains(Vector3(1.0, 0.0, 0.0)), true);   // vertex
        shouldEqual(simplex3.contains(Vector3(0.5, 0.5, 0.0)), true);   // edge
        shouldEqual(simplex3.contains(Vector3(0.25, 0.25, 0.25)), true);    // surface
        shouldEqual(simplex3.contains(Vector3(1.0, 1.0, 1.0)), false);  // outer

        Facet3 facet3(simplex3, 0);
        shouldEqual(facet3.size(), 3);
        shouldEqual(facet3.nSurface(), 0.5);
    }

};

struct SimplexTestSuite : public vigra::test_suite
{
    SimplexTestSuite()
        : vigra::test_suite("SimplexTestSuite")
    {
        add(testCase(&SimplexTest::testBasics));
    }
};

int main(int argc, char** argv)
{
    SimplexTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cerr << test.report() << std::endl;

    return failed != 0;
}

