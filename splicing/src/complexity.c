
#include "splicing.h"
#include "splicing_error.h"

#include <float.h>

int splicing_gene_complexity(const splicing_gff_t *gff, size_t gene,
			     int readLength, int overHang,
			     splicing_complexity_t type,
			     splicing_norm_t norm, int paired,
			     int fast_assignment_matrix,
			     const splicing_vector_t *fragmentProb,
			     int fragmentStart, double normalMean, 
			     double normalVar, double numDevs,
			     double *complexity) {
  
  splicing_matrix_t assignment_matrix;

  SPLICING_CHECK(splicing_matrix_init(&assignment_matrix, 0, 0));
  SPLICING_FINALLY(splicing_matrix_destroy, &assignment_matrix);

  if (!paired) {
    SPLICING_CHECK(splicing_assignment_matrix(gff, gene, readLength, 
					      overHang, &assignment_matrix));
  } else {
    if (!fast_assignment_matrix) {
      SPLICING_CHECK(splicing_paired_assignment_matrix(gff, gene, readLength, 
						       overHang, 
						       fragmentProb, 
						       fragmentStart,
						       normalMean, normalVar,
						       numDevs, 
						       &assignment_matrix));
    } else {
      SPLICING_CHECK(splicing_paired_assignment_matrix_fast(gff, gene,
						       readLength, 
						       overHang, 
						       fragmentProb, 
						       fragmentStart,
						       normalMean, 
						       normalVar, numDevs,
						       &assignment_matrix));
    }
  }

  switch (type) {
  case SPLICING_COMPLEXITY_RELATIVE:
    switch (norm) {
      splicing_vector_t values;
      int i, n;

    case SPLICING_NORM_2:

      SPLICING_CHECK(splicing_vector_init(&values, 0));
      SPLICING_FINALLY(splicing_vector_destroy, &values);
      SPLICING_CHECK(splicing_dgesdd(&assignment_matrix, &values));
      n=splicing_vector_size(&values);
      for (i=n-1; i>=0 && VECTOR(values)[i] < 1e-14; i--) ;
      *complexity = VECTOR(values)[0] / VECTOR(values)[i];
      splicing_vector_destroy(&values);
      SPLICING_FINALLY_CLEAN(1);
      break;

    case SPLICING_NORM_1:

      SPLICING_ERROR("One norm not implemented", SPLICING_UNIMPLEMENTED);
      break;

    case SPLICING_NORM_INFINITY:

      SPLICING_ERROR("Infinity norm not implemented", SPLICING_UNIMPLEMENTED);
      break;

    }
    break;
  case SPLICING_COMPLEXITY_ABSOLUTE:
    SPLICING_ERROR("Absolute complexity not implemented", 
		   SPLICING_UNIMPLEMENTED);
    break;
  }

  splicing_matrix_destroy(&assignment_matrix);
  SPLICING_FINALLY_CLEAN(1);

  return 0;
}

int splicing_cr_matrix(const splicing_gff_t *gff, size_t gene,
		       const splicing_matrix_t *ass_mat, int read_len,
		       int overhang, const splicing_vector_t *expr,
		       int paired, int fast_ass_mat,
		       const splicing_vector_t *fragment_prob,
		       int fragment_start, double normal_mean,
		       double normal_var, double num_devs,
		       splicing_matrix_t *result) {

  splicing_matrix_t r_ass_mat, *my_ass_mat = (splicing_matrix_t *) ass_mat;
  size_t no_iso, no_classes, i, j, k1, k2;
  splicing_vector_t eff_len, X;
  double Y = 0.0;
  double eff_len_no_iso;

  /* Query assignment matrix if not supplied */

  if (! my_ass_mat) {
    my_ass_mat = &r_ass_mat;
    SPLICING_CHECK(splicing_matrix_init(my_ass_mat, 0, 0));
    SPLICING_FINALLY(splicing_matrix_destroy, my_ass_mat);

    if (!paired) {
      SPLICING_CHECK(splicing_assignment_matrix(gff, gene, read_len,
						overhang, my_ass_mat));
    } else if (!fast_ass_mat) {
      SPLICING_CHECK(splicing_paired_assignment_matrix(gff, gene, read_len,
	       overhang, fragment_prob, fragment_start, normal_mean,
	       normal_var, num_devs, my_ass_mat));
    } else {
      SPLICING_CHECK(splicing_paired_assignment_matrix_fast(gff, gene,
	       read_len, overhang, fragment_prob, fragment_start,
	       normal_mean, normal_var, num_devs, my_ass_mat));
    }
  }

  /* Normalize assignment matrix */

  splicing_gff_noiso_one(gff, gene, &no_iso);
  no_classes = splicing_matrix_ncol(my_ass_mat);

  SPLICING_VECTOR_INIT_FINALLY(&eff_len, no_iso);
  for (i = 0; i < no_iso; i++) {
    for (j = 0; j < no_classes; j++) {
      VECTOR(eff_len)[i] += MATRIX(*my_ass_mat, i, j);
    }
  }
  for (i = 0; i < no_iso; i++) {
    for (j = 0; j < no_classes; j++) {
      MATRIX(*my_ass_mat, i, j) /= VECTOR(eff_len)[i];
    }
  }

  SPLICING_VECTOR_INIT_FINALLY(&X, no_classes);
  for (i = 0; i < no_iso; i++) {
    double tmp = VECTOR(eff_len)[i] * VECTOR(*expr)[i];
    Y += tmp;
    for (j = 0; j < no_classes; j++) {
      VECTOR(X)[j] += MATRIX(*my_ass_mat, i, j) * tmp;
    }
  }

  SPLICING_CHECK(splicing_matrix_resize(result, no_iso-1, no_iso-1));

  eff_len_no_iso = VECTOR(eff_len)[no_iso-1];
  for (k1 = 0; k1 < no_iso - 1; k1++) {
    double eff_len_k1 = VECTOR(eff_len)[k1];
    for (k2 = k1; k2 < no_iso - 1; k2++) {
      double C = 0.0;
      double eff_len_k2 = VECTOR(eff_len)[k2];
      for (j = 0; j < no_classes; j++) {
	double Xj = VECTOR(X)[j];
	if (VECTOR(X)[j] > DBL_EPSILON) {
	  double t1 = ((eff_len_k1 - eff_len_no_iso) *
		       (eff_len_k2 - eff_len_no_iso)) / (Y * Y);
	  double t2 = ((MATRIX(*my_ass_mat, k1, j) * eff_len_k1 -
			MATRIX(*my_ass_mat, no_iso-1, j) * eff_len_no_iso) *
		       (MATRIX(*my_ass_mat, k2, j) * eff_len_k2 -
			MATRIX(*my_ass_mat, no_iso-1, j) * eff_len_no_iso))
	               / (Xj * Xj);
	  C += Xj / Y * (t1 - t2);
	}
      }
      /* TODO: handle infinity */
      MATRIX(*result, k1, k2) = MATRIX(*result, k2, k1) = -C;
    }
  }

  splicing_vector_destroy(&X);
  splicing_vector_destroy(&eff_len);
  SPLICING_FINALLY_CLEAN(2);
  if (! ass_mat) {
    splicing_matrix_destroy(my_ass_mat);
    SPLICING_FINALLY_CLEAN(1);
  }

  return 0;
}
