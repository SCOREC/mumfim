#ifndef MUMFIM_SRC_MUMFIM_MICROSCALE_TENSORUTILITIES_H
#define MUMFIM_SRC_MUMFIM_MICROSCALE_TENSORUTILITIES_H
#include <Kokkos_Core.hpp>

namespace mumfim
{

  /**
   * @brief Compute the determinant of a 2x2 or 3x3 matrix
   * @param[in] F the matrix
   * @return the determinant
   * @note this function is KOKKOS_INLINE and can be used on the device
   */
  template <int dim, typename ViewT>
  KOKKOS_INLINE_FUNCTION constexpr auto Determinant(ViewT F)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 2,
                          double>
  {
    static_assert(dim == 2 || dim == 3, "only 2d and 3d implemented");
    assert(F.extent(0) == dim && F.extent(1) == dim);
    if constexpr (dim == 2)
    {
      return (F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0));
    }
    else if constexpr (dim == 3)
    {
      return (F(0, 0) * F(1, 1) * F(2, 2) + F(0, 1) * F(1, 2) * F(2, 0) +
              F(0, 2) * F(1, 0) * F(2, 1) -
              (F(0, 2) * F(1, 1) * F(2, 0) + F(0, 1) * F(1, 0) * F(2, 2) +
               F(0, 0) * F(1, 2) * F(2, 1)));
    }
  }

  /**
   * @brief Compute the inverse of a 2x2 or 3x3 matrix
   * @param[in] F the matrix
   * @param[out] F_inv the inverse of F
   * @return the determinant of F
   * @note this function is KOKKOS_INLINE and can be used on the device
   */
  template <int dim, typename ViewT>
  KOKKOS_INLINE_FUNCTION constexpr auto Invert(ViewT F, ViewT invF)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 2,
                          void>
  {
    static_assert(dim == 2 || dim == 3, "only 2d and 3d implemented");
    assert(F.extent(0) == dim && F.extent(1) == dim);
    if constexpr (dim == 2)
    {
      double detF = Determinant<2>(F);
      KOKKOS_ASSERT(detF > 0.0);
      invF(0, 0) = F(1, 1) / detF;
      invF(0, 1) = -F(0, 1) / detF;
      invF(1, 0) = -F(1, 0) / detF;
      invF(1, 1) = F(0, 0) / detF;
    }
    else if constexpr (dim == 3)
    {
      double detF = Determinant<3>(F);
      KOKKOS_ASSERT(detF > 0.0);
      invF(0, 0) = (F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1)) / detF;
      invF(0, 1) = (F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2)) / detF;
      invF(0, 2) = (F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1)) / detF;
      invF(1, 0) = (F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2)) / detF;
      invF(1, 1) = (F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0)) / detF;
      invF(1, 2) = (F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2)) / detF;
      invF(2, 0) = (F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0)) / detF;
      invF(2, 1) = (F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1)) / detF;
      invF(2, 2) = (F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0)) / detF;
    }
  }

  template <typename ExeSpace>
  void MandelToVoigt(
      Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space> matrix)
  {
    Kokkos::parallel_for(
        "mandel to voigt",
        Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>>(
            {0, 0, 0}, {matrix.extent(0), matrix.extent(1), matrix.extent(2)}),
        KOKKOS_LAMBDA(int i, int j, int k) {
          auto sr2 = sqrt(2);
          // lower right
          if (j > 2 && k > 2)
          {
            matrix(i, j, k) = matrix(i, j, k) / 2;
          }
          // upper right
          else if (j < 3 && k > 2)
          {
            matrix(i, j, k) = matrix(i, j, k) / sr2;
          }
          // lower left
          else if (j > 2)
          {
            matrix(i, j, k) = matrix(i, j, k) / sr2;
          }
        });
  }

  template <typename ExeSpace>
  void VoigtToMandel(
      Kokkos::View<Scalar * [6][6], typename ExeSpace::memory_space> matrix)
  {
    Kokkos::parallel_for(
        "voigt to mandel",
        Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<3>,
                              Kokkos::IndexType<size_t>>(
            {0ul, 0ul, 0ul},
            {matrix.extent(0), matrix.extent(1), matrix.extent(2)}),
        KOKKOS_LAMBDA(int i, int j, int k) {
          auto sr2 = sqrt(2);
          // lower right
          if (j > 2 && k > 2)
          {
            matrix(i, j, k) = 2 * matrix(i, j, k);
          }
          // upper right
          else if (j < 3 && k > 2)
          {
            matrix(i, j, k) = sr2 * matrix(i, j, k);
          }
          // lower left
          else if (j > 2)
          {
            matrix(i, j, k) = sr2 * matrix(i, j, k);
          }
        });
  }

  /**
   * @brief Convert a stress tensor from Voigt to matrix form
   * @param[in] stress_voigt the stress tensor in Voigt form
   * @param[out] stress_matrix the stress tensor in matrix form
   * @note this function is KOKKOS_INLINE and can be used on the device
   */
  template <typename ViewT, typename ViewT2>
  KOKKOS_INLINE_FUNCTION constexpr auto VoigtToMatrix(ViewT stress_voigt,
                                                      ViewT2 stress_matrix)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 1 &&
                              Kokkos::is_view<ViewT2>::value &&
                              ViewT2::rank == 2,
                          void>
  {
    KOKKOS_ASSERT(stress_voigt.extent(0) == 6);
    KOKKOS_ASSERT(stress_matrix.extent(0) == 3 && stress_matrix.extent(1) == 3);
    stress_matrix(0, 0) = stress_voigt(0);
    stress_matrix(1, 1) = stress_voigt(1);
    stress_matrix(2, 2) = stress_voigt(2);
    stress_matrix(2, 1) = stress_voigt(3);
    stress_matrix(1, 2) = stress_voigt(3);
    stress_matrix(0, 2) = stress_voigt(4);
    stress_matrix(2, 0) = stress_voigt(4);
    stress_matrix(0, 1) = stress_voigt(5);
    stress_matrix(1, 0) = stress_voigt(5);
  }

  /**
   * @brief Convert a stress tensor from matrix to Voigt form
   * @param[in] stress_matrix the stress tensor in matrix form
   * @param[out] stress_voigt the stress tensor in Voigt form
   * @note this function is KOKKOS_INLINE and can be used on the device
   */
  template <typename ViewT, typename ViewT2>
  KOKKOS_INLINE_FUNCTION constexpr auto MatrixToVoigt(ViewT stress_matrix,
                                                      ViewT2 stress_voigt)
      -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 2 &&
                              Kokkos::is_view<ViewT2>::value &&
                              ViewT2::rank == 1,
                          void>
  {
    stress_voigt(0) = stress_matrix(0, 0);
    stress_voigt(1) = stress_matrix(1, 1);
    stress_voigt(2) = stress_matrix(2, 2);
    stress_voigt(3) = stress_matrix(2, 1);
    stress_voigt(4) = stress_matrix(0, 2);
    stress_voigt(5) = stress_matrix(1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  size_t TensorIndex2VoigtIndex(size_t i, size_t j)
  {
    // 11, 22, 33, 23, 13, 12
    if (i == j)
    {
      return i;
    }
    else if ((i == 1 && j == 2) || (i == 2 && j == 1))
    {
      return 3;
    }
    else if ((i == 0 && j == 2) || (i == 2 && j == 0))
    {
      return 4;
    }
    else if ((i == 0 && j == 1) || (i == 1 && j == 0))
    {
      return 5;
    }
    else
    {
      KOKKOS_ASSERT(i < 3 && j < 3);
    }
    return -1;
  }
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_MICROSCALE_TENSORUTILITIES_H
