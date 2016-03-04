#ifndef VIGRA_SIMPLEX_HXX
#define VIGRA_SIMPLEX_HXX

#include <cmath>
#include <algorithm>
#include "config.hxx"
#include "error.hxx"
#include "tinyvector.hxx"
#include "array_vector.hxx"
#include "linear_algebra.hxx"

namespace vigra {

template <class POINT> class Simplex;

template <class POINT>
class Facet
: protected ArrayVector<POINT>
{
  public:
    enum Dimension {dimension = POINT::static_size - 1};
    typedef ArrayVector<POINT> Base;
    typedef POINT                                 Point;
    typedef typename NumericTraits<Point>::RealPromote RealPoint;
    typedef typename Base::value_type             value_type;
    typedef typename Base::reference              reference;
    typedef typename Base::const_reference        const_reference;
    typedef typename Base::pointer                pointer;
    typedef typename Base::const_pointer          const_pointer;
    typedef typename Base::iterator               iterator;
    typedef typename Base::const_iterator         const_iterator;
    typedef typename Base::reverse_iterator       reverse_iterator;
    typedef typename Base::const_reverse_iterator const_reverse_iterator;
    typedef typename Base::size_type              size_type;
    typedef typename Base::difference_type        difference_type;
    typedef typename Point::value_type            coordinate_type;
    typedef typename RealPoint::value_type        real_type;

    using Base::size;
    using Base::empty;
    using Base::begin;
    using Base::end;
    using Base::cbegin;
    using Base::cend;
    using Base::rbegin;
    using Base::rend;
    using Base::crbegin;
    using Base::crend;

    Facet()
    : Base(dimension + 1, Point())
    , normal_(Point())
    {};
    
    Facet(const Simplex<Point> & simplex, difference_type index)
    : Base(dimension + 1, Point())
    , normal_(Point())
    {
        vigra_precondition(index < simplex.size(),
                "Facet::Facet() index must be smaller than simplex size");
        // copy set of facet vertices
        for (int i = 0, j=0; i+1 < simplex.size(); i++, j++)
        {
            if (i == index)
            {
                j++;
            }
            (*this)[i] = simplex[j];
        }
        // get the normal vector
        MultiArray<2, real_type> mat(dimension+1, dimension+1);
        {
            Point vec = (*this)[0] - simplex[index];
            std::copy(vec.begin(),
                      vec.end(),
                      mat.template bind<0>(dimension).begin());
        }
        for (int i = 0; i < size() - 1; i++)
        {
            Point vec = (*this)[i+1] - (*this)[0];
            std::copy(vec.begin(), vec.end(), mat.template bind<0>(i).begin());
        }
        linalg::inverse(mat, mat);
        std::copy(mat.template bind<1>(dimension).begin(),
                  mat.template bind<1>(dimension).end(),
                  normal_.begin());
    }

    const RealPoint& normal() const
    {
        return normal_;
    }

    double nSurface() const
    {
        MultiArray<2, coordinate_type> mat(Shape2(dimension + 1, dimension + 1));
        double fac = 1.0;
        for (int i = 0; i < dimension; i++) {
            fac *= i+1;
            Point vec = (*this)[i+1] - (*this)[0];
            std::copy(vec.begin(), vec.end(), mat.template bind<0>(i).begin());
        }
        std::copy(normal_.begin(),
                  normal_.end(),
                  mat.template bind<0>(dimension).begin());
        return abs(linalg::determinant(mat) / fac);
    }

  private:
    RealPoint normal_;
};

template <class POINT>
class Simplex
: protected ArrayVector<POINT>
{
  public:
    enum Dimension {dimension = POINT::static_size};
    typedef ArrayVector<POINT> Base;
    typedef POINT                                 Point;
    typedef typename Base::value_type             value_type;
    typedef typename Base::reference              reference;
    typedef typename Base::const_reference        const_reference;
    typedef typename Base::pointer                pointer;
    typedef typename Base::const_pointer          const_pointer;
    typedef typename Base::iterator               iterator;
    typedef typename Base::const_iterator         const_iterator;
    typedef typename Base::reverse_iterator       reverse_iterator;
    typedef typename Base::const_reverse_iterator const_reverse_iterator;
    typedef typename Base::size_type              size_type;
    typedef typename Base::difference_type        difference_type;
    typedef typename POINT::value_type            coordinate_type;

    using Base::size;
    using Base::empty;
    using Base::begin;
    using Base::end;
    using Base::cbegin;
    using Base::cend;
    using Base::rbegin;
    using Base::rend;
    using Base::crbegin;
    using Base::crend;
    using Base::operator[];

    Simplex()
    : Base(dimension, Point())
    {}

    template <class InputIterator>
    Simplex(InputIterator b, InputIterator e)
    : Base(b, e)
    {
        vigra_precondition(e - b == dimension + 1,
                "Simplex::Simplex() Number of vectors must match dimensions");
    }

    double nVolume() const
    {
        MultiArray<2, coordinate_type> mat(Shape2(dimension, dimension));
        double fac = 1.0;
        for (int n = 1; n < dimension+1; n++)
        {
            fac *= n;
            Point vec((*this)[n] - (*this)[0]);
            std::copy(vec.begin(), vec.end(), mat.template bind<0>(n-1).begin());
        }
        return abs(linalg::determinant(mat) / fac);
    }

    double nSurface() const
    {
        double ret = 0;
        for (int i = 0; i < size(); i++)
        {
            Facet<Point> facet((*this), i);
            ret += facet.nSurface();
        }
        return ret;
    }

    bool contains(const_reference point) const
    {
        MultiArray<2, coordinate_type> jp_mat(Shape2(dimension, dimension));
        MultiArray<2, coordinate_type> jj_mat(Shape2(dimension, dimension));
        for (int j = 0; j < dimension + 1; j++)
        {
            for (int i = 0, ii = 0; i < dimension; i++, ii++)
            {
                if (i == j)
                {
                    ii++;
                }
                Point u((*this)[ii] - point);
                std::copy(u.begin(),
                          u.end(), 
                          jp_mat.template bind<0>(i).begin());
                Point v((*this)[ii] - (*this)[j]);
                std::copy(v.begin(),
                          v.end(),
                          jj_mat.template bind<0>(i).begin());
            }
            coordinate_type jj_det = linalg::determinant(jj_mat);
            coordinate_type jp_det = linalg::determinant(jp_mat);
            if (((jj_det > 0) xor (jp_det > 0)) and (jp_det != 0))
            {
                return false;
            }
        }
        return true;
    }
};

} // namespace vigra

#endif
