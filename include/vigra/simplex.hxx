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

template <unsigned int N, class T> class SimplexView;
template <unsigned int N, class T> class Simplex;

template <unsigned int N, class T>
class FacetView
: protected ArrayVector<TinyVectorView<T, N> >
{
  public:
    enum Dimension {dimension = N};
    typedef ArrayVector<TinyVectorView<T, N> >         base_type;
    typedef TinyVectorView<T, N>                       point_type;
    typedef typename base_type::value_type             value_type;
    typedef typename base_type::reference              reference;
    typedef typename base_type::const_reference        const_reference;
    typedef typename base_type::pointer                pointer;
    typedef typename base_type::const_pointer          const_pointer;
    typedef typename base_type::iterator               iterator;
    typedef typename base_type::const_iterator         const_iterator;
    typedef typename base_type::reverse_iterator       reverse_iterator;
    typedef typename base_type::const_reverse_iterator const_reverse_iterator;
    typedef typename base_type::size_type              size_type;
    typedef typename base_type::difference_type        difference_type;
    typedef typename point_type::value_type            coordinate_type;
    typedef typename NumericTraits<T>::RealPromote     real_type;
    typedef TinyVector<real_type, N>                   real_point_type;

    using base_type::size;
    using base_type::empty;
    using base_type::begin;
    using base_type::end;
    using base_type::cbegin;
    using base_type::cend;
    using base_type::rbegin;
    using base_type::rend;
    using base_type::crbegin;
    using base_type::crend;

    FacetView()
    : base_type()
    , normal_(point_type())
    {};
    
    FacetView(const SimplexView<N, T> & simplex, difference_type index)
    : base_type()
    , normal_(real_point_type())
    {
        vigra_precondition(index < simplex.size(),
                "Facet::Facet() Index must be smaller than simplex size");
        // copy set of facet vertices
        for (int i = 0, j=0; i+1 < simplex.size(); i++, j++)
        {
            if (i == index)
            {
                j++;
            }
            this->push_back(simplex[j]);
        }
        calculateNormal(simplex[index]);
    }

    template <class InputIterator>
    FacetView(InputIterator b, InputIterator e)
    : base_type(b, e - 1)
    , normal_(real_point_type())
    {
        vigra_precondition(e - b == dimension + 1,
                           "Facet::Facet() Number of points must be dim + 1");
        calculateNormal(*(e-1));
    }

    void calculateNormal(const point_type & inner)
    {
        // get the normal vector
        MultiArray<2, real_type> mat(dimension, dimension);
        mat.template bind<0>(dimension-1) = (*this)[0] - inner;
        for (int i = 0; i < size() - 1; i++)
        {
            mat.template bind<0>(i) = (*this)[i+1] - (*this)[0];
        }
        linalg::inverse(mat, mat);
        std::copy(mat.template bind<1>(dimension-1).begin(),
                  mat.template bind<1>(dimension-1).end(),
                  normal_.begin());
        normal_ /= norm(normal_);
    }

    virtual const real_point_type& normal() const
    {
        return normal_;
    }

    virtual double nSurface() const
    {
        MultiArray<2, real_type> mat(Shape2(dimension, dimension));
        double fac = 1.0;
        for (int i = 0; i < dimension - 1; i++)
        {
            fac *= i+1;
            mat.template bind<0>(i) = (*this)[i+1] - (*this)[0];
        }
        mat.template bind<0>(dimension - 1) = normal_;
        return abs(linalg::determinant(mat) / fac);
    }

  protected:
    real_point_type normal_;
};

template <unsigned int N, class T>
class SimplexView
: protected ArrayVector<TinyVectorView<T, N> >
{
  friend class Simplex<N, T>;

  public:
    enum Dimension {dimension = N};
    typedef ArrayVector<TinyVectorView<T, N> >         base_type;
    typedef typename base_type::value_type             point_type;
    typedef typename base_type::value_type             value_type;
    typedef typename base_type::reference              reference;
    typedef typename base_type::const_reference        const_reference;
    typedef typename base_type::pointer                pointer;
    typedef typename base_type::const_pointer          const_pointer;
    typedef typename base_type::iterator               iterator;
    typedef typename base_type::const_iterator         const_iterator;
    typedef typename base_type::reverse_iterator       reverse_iterator;
    typedef typename base_type::const_reverse_iterator const_reverse_iterator;
    typedef typename base_type::size_type              size_type;
    typedef typename base_type::difference_type        difference_type;
    typedef typename point_type::value_type            coordinate_type;
    typedef FacetView<N, T>                            facet_type;

    using base_type::size;
    using base_type::empty;
    using base_type::begin;
    using base_type::end;
    using base_type::cbegin;
    using base_type::cend;
    using base_type::rbegin;
    using base_type::rend;
    using base_type::crbegin;
    using base_type::crend;
    using base_type::operator[];

  protected:
    using base_type::push_back;

  public:

    SimplexView()
    : base_type()
    , eps_(std::numeric_limits<T>::epsilon() * 2)
    {}

    template <class InputIterator>
    SimplexView(InputIterator b, InputIterator e)
    : base_type(b, e)
    , eps_(std::numeric_limits<T>::epsilon() * 2)
    {
        vigra_precondition(e - b == dimension + 1,
                "Simplex::Simplex() Number of vectors must be dim + 1");
    }

    virtual double nVolume() const
    {
        MultiArray<2, coordinate_type> mat(Shape2(dimension, dimension));
        double fac = 1.0;
        for (int n = 1; n < dimension+1; n++)
        {
            fac *= n;
            mat.template bind<0>(n-1) = (*this)[n] - (*this)[0];
        }
        return abs(linalg::determinant(mat) / fac);
    }

    virtual double nSurface() const
    {
        double ret = 0;
        for (int n = 0; n < dimension+1; n++)
        {
            facet_type facet(*this, n);
            ret += facet.nSurface();
        }
        return ret;
    }

    virtual bool contains(const_reference point) const
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
                jp_mat.template bind<0>(i) = (*this)[ii] - point;
                jj_mat.template bind<0>(i) = (*this)[ii] - (*this)[j];
            }
            coordinate_type jj_det = linalg::determinant(jj_mat);
            coordinate_type jp_det = linalg::determinant(jp_mat);
            if (((jj_det > 0) xor (jp_det > 0)) and abs(jp_det) > eps_)
            {
                return false;
            }
        }
        return true;
    }

    const T eps_;
};

template <unsigned int N, class T>
class Simplex
: public SimplexView<N, T>
{
  public:
    enum Dimension {dimension = N};
    typedef SimplexView<N, T>                          base_type;
    typedef typename base_type::value_type             point_type;
    typedef typename base_type::value_type             value_type;
    typedef typename base_type::reference              reference;
    typedef typename base_type::const_reference        const_reference;
    typedef typename base_type::pointer                pointer;
    typedef typename base_type::const_pointer          const_pointer;
    typedef typename base_type::iterator               iterator;
    typedef typename base_type::const_iterator         const_iterator;
    typedef typename base_type::reverse_iterator       reverse_iterator;
    typedef typename base_type::const_reverse_iterator const_reverse_iterator;
    typedef typename base_type::size_type              size_type;
    typedef typename base_type::difference_type        difference_type;
    typedef typename point_type::value_type            coordinate_type;
    typedef FacetView<N, T>                            facet_type;

    using base_type::size;
    using base_type::empty;
    using base_type::begin;
    using base_type::end;
    using base_type::cbegin;
    using base_type::cend;
    using base_type::rbegin;
    using base_type::rend;
    using base_type::crbegin;
    using base_type::crend;
    using base_type::operator[];

    Simplex()
    : base_type()
    {}

    Simplex(const MultiArrayView<2, T> & mat)
    : base_type()
    {
        vigra_precondition(mat.shape(0) == dimension + 1,
                           "Simplex::Simplex(): Shape of input does not match");
        vigra_precondition(mat.shape(1) == dimension,
                           "Simplex::Simplex(): Shape of input does not macht");
        std::copy(mat.begin(), mat.end(), data_.begin());
        for (int n = 0; n < N+1; n++)
        {
            this->push_back(point_type(&(data_[n*N])));
        }
    }

    template <class InputIterator>
    Simplex(InputIterator b, InputIterator e)
    : base_type()
    {
        vigra_precondition(e - b == dimension + 1,
                           "Simplex::Simplex(): Number of vertices must be dim+1");
        InputIterator it = b;
        for (int n = 0; n < dimension + 1; n++, it++)
        {
            std::copy(it->begin(), it->end(), &(data_[n*N]));
            this->push_back(point_type(&(data_[n*N])));
        }
    }

  private:
    std::array<T, N*(N+1)> data_;
};

} // namespace vigra

#endif
