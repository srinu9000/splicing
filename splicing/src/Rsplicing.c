
#include <R.h>
#include <Rdefines.h>

#include "splicing_random.h"
#include "splicing_error.h"
#include "Rsplicing.h"

void R_splicing_error_handler(const char *reason, const char *file, 
			      int line, int splicing_errno) {
  SPLICING_FINALLY_FREE();
  error("At %s:%i : %s, %s", file, line, reason,
	splicing_strerror(splicing_errno));
}

void R_splicing_warning_handler(const char *reason, const char *file,
				int line, int splicing_errno) {
  warning("At %s:%i : %s, %s", file, line, reason,
	  splicing_strerror(splicing_errno));
}

int R_splicing_begin() {
  splicing_set_error_handler(R_splicing_error_handler);
  splicing_set_warning_handler(R_splicing_warning_handler);
  GetRNGstate();
  return 0;
}

int R_splicing_end() {
  PutRNGstate();
  return 0;
}

SEXP R_splicing_finalizer() {
  SPLICING_FINALLY_FREE();
  return R_NilValue;
}

SEXP R_splicing_getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

SEXP R_splicing_vector_int_to_SEXP(const splicing_vector_int_t *v) {
  SEXP result;
  size_t n=splicing_vector_int_size(v);
  
  PROTECT(result=NEW_INTEGER(n));
  memcpy(INTEGER(result), VECTOR(*v), n*sizeof(int));
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_vector_to_SEXP(const splicing_vector_t *v) {
  SEXP result;
  size_t n=splicing_vector_size(v);
  
  PROTECT(result=NEW_NUMERIC(n));
  memcpy(REAL(result), VECTOR(*v), n*sizeof(double));
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_matrix_to_SEXP(const splicing_matrix_t *m) {
  int r = splicing_matrix_nrow(m);
  int c = splicing_matrix_ncol(m);
  int n = r * c;
  SEXP result, dim;
  
  PROTECT(result = NEW_NUMERIC(n));
  PROTECT(dim =  NEW_INTEGER(2));
  INTEGER(dim)[0] = r;
  INTEGER(dim)[1] = c;
  SET_DIM(result, dim);

  memcpy(REAL(result), &MATRIX(*m, 0, 0), n * sizeof(double));
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_matrix_int_to_SEXP(const splicing_matrix_int_t *m) {
  int r = splicing_matrix_int_nrow(m);
  int c = splicing_matrix_int_ncol(m);
  int n = r * c;
  SEXP result, dim;
  
  PROTECT(result = NEW_INTEGER(n));
  PROTECT(dim =  NEW_INTEGER(2));
  INTEGER(dim)[0] = r;
  INTEGER(dim)[1] = c;
  SET_DIM(result, dim);

  memcpy(INTEGER(result), &MATRIX(*m, 0, 0), n * sizeof(int));
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_vector_int_index_to_SEXP(const splicing_vector_int_t *v,
					 const splicing_vector_int_t *idx) {
  int i, n=splicing_vector_int_size(idx);
  int vn=splicing_vector_int_size(v);
  SEXP result;
  
  PROTECT(result=NEW_LIST(n));

  for (i=0; i<n; i++) {
    int j=VECTOR(*idx)[i];
    int x, k= i<n-1 ? VECTOR(*idx)[i+1] : vn;
    SEXP tmp;
    int *itmp;
    PROTECT(tmp=NEW_INTEGER(k-j));
    itmp=INTEGER(tmp);
    for (x=0; j<k; j++, x++) {
      itmp[x] = VECTOR(*v)[j];
    }
    SET_VECTOR_ELT(result, i, tmp);
    UNPROTECT(1);
  }
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_miso_rundata_to_SEXP(const splicing_miso_rundata_t *data) {
  SEXP result, names;
  
  PROTECT(result=NEW_LIST(8));
  SET_VECTOR_ELT(result, 0, ScalarInteger(data->noIso));
  SET_VECTOR_ELT(result, 1, ScalarInteger(data->noIters));
  SET_VECTOR_ELT(result, 2, ScalarInteger(data->noBurnIn));
  SET_VECTOR_ELT(result, 3, ScalarInteger(data->noLag));
  SET_VECTOR_ELT(result, 4, ScalarInteger(data->noAccepted));
  SET_VECTOR_ELT(result, 5, ScalarInteger(data->noRejected));
  SET_VECTOR_ELT(result, 6, ScalarInteger(data->noChains));
  SET_VECTOR_ELT(result, 7, ScalarInteger(data->noSamples));
  
  PROTECT(names=NEW_STRING(8));
  SET_STRING_ELT(names, 0, mkChar("noIso"));
  SET_STRING_ELT(names, 1, mkChar("noIters"));
  SET_STRING_ELT(names, 2, mkChar("noBurnIn"));
  SET_STRING_ELT(names, 3, mkChar("noLag"));
  SET_STRING_ELT(names, 4, mkChar("noAccepted"));
  SET_STRING_ELT(names, 5, mkChar("noRejected"));
  SET_STRING_ELT(names, 6, mkChar("noChains"));
  SET_STRING_ELT(names, 7, mkChar("noSamples"));
  
  SET_NAMES(result, names);
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_reads_to_SEXP(const splicing_reads_t *reads) {
  SEXP result, names, class;
  
  PROTECT(result=NEW_LIST(18));
  SET_VECTOR_ELT(result, 0, R_splicing_strvector_to_SEXP(&reads->chrname));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&reads->chrlen));
  SET_VECTOR_ELT(result, 2, R_splicing_vector_int_to_SEXP(&reads->chr));
  SET_VECTOR_ELT(result, 3, R_splicing_strvector_to_SEXP(&reads->qname));
  SET_VECTOR_ELT(result, 4, R_splicing_strvector_to_SEXP(&reads->cigar));
  SET_VECTOR_ELT(result, 5, R_splicing_vector_int_to_SEXP(&reads->position));
  SET_VECTOR_ELT(result, 6, R_splicing_vector_int_to_SEXP(&reads->flags));
  SET_VECTOR_ELT(result, 7, R_splicing_vector_int_to_SEXP(&reads->pairpos));
  SET_VECTOR_ELT(result, 8, ScalarInteger(reads->noPairs));
  SET_VECTOR_ELT(result, 9, ScalarInteger(reads->noSingles));
  SET_VECTOR_ELT(result, 10, ScalarLogical(reads->paired));
  SET_VECTOR_ELT(result, 11, R_splicing_vector_int_to_SEXP(&reads->mapq));
  SET_VECTOR_ELT(result, 12, R_splicing_vector_int_to_SEXP(&reads->rnext));
  SET_VECTOR_ELT(result, 13, R_splicing_vector_int_to_SEXP(&reads->tlen));
  SET_VECTOR_ELT(result, 14, R_splicing_strvector_to_SEXP(&reads->seq));
  SET_VECTOR_ELT(result, 15, R_splicing_strvector_to_SEXP(&reads->qual));
  SET_VECTOR_ELT(result, 16, R_splicing_vector_int_to_SEXP(&reads->mypair));
  SET_VECTOR_ELT(result, 17, R_splicing_strvector_to_SEXP(&reads->attributes));

  PROTECT(names=NEW_CHARACTER(18));
  SET_STRING_ELT(names, 0, mkChar("chrname"));
  SET_STRING_ELT(names, 1, mkChar("chrlen"));
  SET_STRING_ELT(names, 2, mkChar("chr"));
  SET_STRING_ELT(names, 3, mkChar("qname"));
  SET_STRING_ELT(names, 4, mkChar("cigar"));
  SET_STRING_ELT(names, 5, mkChar("position"));
  SET_STRING_ELT(names, 6, mkChar("flag"));
  SET_STRING_ELT(names, 7, mkChar("pairpos"));
  SET_STRING_ELT(names, 8, mkChar("noPairs"));
  SET_STRING_ELT(names, 9, mkChar("noSingles"));
  SET_STRING_ELT(names, 10, mkChar("paired"));
  SET_STRING_ELT(names, 11, mkChar("mapq"));
  SET_STRING_ELT(names, 12, mkChar("rnext"));
  SET_STRING_ELT(names, 13, mkChar("tlen"));
  SET_STRING_ELT(names, 14, mkChar("seq"));
  SET_STRING_ELT(names, 15, mkChar("qual"));
  SET_STRING_ELT(names, 16, mkChar("mypair"));
  SET_STRING_ELT(names, 17, mkChar("attributes"));
  SET_NAMES(result, names);

  PROTECT(class=NEW_CHARACTER(1));
  SET_STRING_ELT(class, 0, mkChar("splicingSAM"));
  SET_CLASS(result, class);

  UNPROTECT(3);
  return result;
}

int R_splicing_SEXP_to_reads(SEXP preads, splicing_reads_t *reads) {

  SEXP chrname, chrlen, chr, qname, cigar, position, flag, pairpos,
    noPairs, noSingles, paired, mapq, rnext, tlen, seq, qual, mypair, 
    attributes;

  chrname = R_splicing_getListElement(preads, "chrname");
  chrlen = R_splicing_getListElement(preads, "chrlen");
  chr = R_splicing_getListElement(preads, "chr");
  qname = R_splicing_getListElement(preads, "qname");
  cigar = R_splicing_getListElement(preads, "cigar");
  position = R_splicing_getListElement(preads, "position");
  flag = R_splicing_getListElement(preads, "flag");
  pairpos = R_splicing_getListElement(preads, "pairpos");
  noPairs = R_splicing_getListElement(preads, "noPairs");
  noSingles = R_splicing_getListElement(preads, "noSingles");
  paired = R_splicing_getListElement(preads, "paired");
  mapq = R_splicing_getListElement(preads, "mapq");
  rnext = R_splicing_getListElement(preads, "rnext");
  tlen = R_splicing_getListElement(preads, "tlen");
  seq = R_splicing_getListElement(preads, "seq");
  qual = R_splicing_getListElement(preads, "qual");
  mypair = R_splicing_getListElement(preads, "mypair");
  attributes = R_splicing_getListElement(preads, "attributes");
  
  R_splicing_SEXP_to_strvector(chrname, &reads->chrname);
  R_splicing_SEXP_to_vector_int(chrlen, &reads->chrlen);
  R_splicing_SEXP_to_vector_int(chr, &reads->chr);
  R_splicing_SEXP_to_strvector(qname, &reads->qname);
  R_splicing_SEXP_to_strvector(cigar, &reads->cigar);
  R_splicing_SEXP_to_vector_int(position, &reads->position);
  R_splicing_SEXP_to_vector_int(pairpos, &reads->pairpos);
  R_splicing_SEXP_to_vector_int(flag, &reads->flags);
  reads->noPairs=INTEGER(noPairs)[0];
  reads->noSingles=INTEGER(noSingles)[0];
  reads->paired=LOGICAL(paired)[0];
  R_splicing_SEXP_to_vector_int(mapq, &reads->mapq);
  R_splicing_SEXP_to_vector_int(rnext, &reads->rnext);
  R_splicing_SEXP_to_vector_int(tlen, &reads->tlen);
  R_splicing_SEXP_to_strvector(seq, &reads->seq);
  R_splicing_SEXP_to_strvector(qual, &reads->qual);
  R_splicing_SEXP_to_vector_int(mypair, &reads->mypair);
  R_splicing_SEXP_to_strvector(attributes, &reads->attributes);
  
  return 0;
}

int R_splicing_SEXP_to_vector_int(SEXP pv, splicing_vector_int_t *v) {
  v->stor_begin = INTEGER(pv);
  v->stor_end = v->end = v->stor_begin + GET_LENGTH(pv);
  return 0;
}

int R_splicing_SEXP_to_vector_int_copy(SEXP pv, splicing_vector_int_t *v) {
  int n=GET_LENGTH(pv);
  splicing_vector_int_init(v, n);
  memcpy(VECTOR(*v), INTEGER(pv), sizeof(int) * n);
  return 0;
}

int R_splicing_SEXP_to_vector(SEXP pv, splicing_vector_t *v) {
  v->stor_begin = REAL(pv);
  v->stor_end = v->end = v->stor_begin + GET_LENGTH(pv);
  return 0;
}

int R_splicing_SEXP_to_vector_copy(SEXP pv, splicing_vector_t *v) {
  int n=GET_LENGTH(pv);
  splicing_vector_init(v, n);
  memcpy(VECTOR(*v), REAL(pv), sizeof(int) * n);
  return 0;
}

int R_splicing_SEXP_to_matrix(SEXP pm, splicing_matrix_t *m) {
  int *dim=INTEGER(GET_DIM(pm));
  R_splicing_SEXP_to_vector(pm, &m->data);
  m->nrow=dim[0];
  m->ncol=dim[1];
  return 0;
}

int R_splicing_SEXP_to_matrix_int(SEXP pm, splicing_matrix_int_t *m) {
  int *dim=INTEGER(GET_DIM(pm));
  R_splicing_SEXP_to_vector_int(pm, &m->data);
  m->nrow=dim[0];
  m->ncol=dim[1];
  return 0;
}

/* We construct a strvector by hand */

int R_splicing_SEXP_to_strvector(SEXP pv, splicing_strvector_t *v) {
  size_t i;
  v->size = v->asize = GET_LENGTH(pv);
  v->free = 0;
  v->table = (char **) R_alloc(v->size, sizeof(char*));
  for (i=0; i<v->size; i++) {
    v->table[i] = (char*) CHAR(STRING_ELT(pv, i));
  }

  return 1;
}

int R_splicing_SEXP_to_strvector_copy(SEXP pv, splicing_strvector_t *v) {
  long int i;
  v->free = 1;
  splicing_strvector_init(v, GET_LENGTH(pv));
  for (i=0; i<v->size; i++) {
    splicing_strvector_set(v, i, CHAR(STRING_ELT(pv, i)));
  }

  return 0;
}

int R_splicing_SEXP_to_gff(SEXP pgff, splicing_gff_t *gff) {

  SEXP seqid, source, genes, transcripts, type, start, end, score, 
    strand, phase, ID, parent, seqid_str, source_str, attributes;

  seqid = R_splicing_getListElement(pgff, "seqid");
  source = R_splicing_getListElement(pgff, "source");
  genes = R_splicing_getListElement(pgff, "gid");
  transcripts = R_splicing_getListElement(pgff, "tid");
  type = R_splicing_getListElement(pgff, "type");
  start = R_splicing_getListElement(pgff, "start");
  end = R_splicing_getListElement(pgff, "end");
  score = R_splicing_getListElement(pgff, "score");
  strand = R_splicing_getListElement(pgff, "strand");
  phase = R_splicing_getListElement(pgff, "phase");
  ID = R_splicing_getListElement(pgff, "ID");
  parent = R_splicing_getListElement(pgff, "parent");
  seqid_str = R_splicing_getListElement(pgff, "seqid_str");
  source_str = R_splicing_getListElement(pgff, "source_str");
  attributes = R_splicing_getListElement(pgff, "attributes");

  R_splicing_SEXP_to_vector_int(seqid, &gff->seqid);
  R_splicing_SEXP_to_vector_int(source, &gff->source);
  R_splicing_SEXP_to_vector_int(genes, &gff->genes);
  R_splicing_SEXP_to_vector_int(transcripts, &gff->transcripts);
  R_splicing_SEXP_to_vector_int(type, &gff->type);
  R_splicing_SEXP_to_vector_int(start, &gff->start);
  R_splicing_SEXP_to_vector_int(end, &gff->end);
  R_splicing_SEXP_to_vector(score, &gff->score);
  R_splicing_SEXP_to_vector_int(strand, &gff->strand);
  R_splicing_SEXP_to_vector_int(phase, &gff->phase);
  R_splicing_SEXP_to_strvector(ID, &gff->ID);
  R_splicing_SEXP_to_vector_int(parent, &gff->parent);
  R_splicing_SEXP_to_strvector(seqid_str, &gff->seqids);
  R_splicing_SEXP_to_strvector(source_str, &gff->sources);

  gff->n = GET_LENGTH(start);
  gff->nogenes = GET_LENGTH(genes);
  gff->notranscripts = GET_LENGTH(transcripts);

  return 0;
}

int R_splicing_SEXP_to_gff_copy(SEXP pgff, splicing_gff_t *gff) {

  SEXP seqid, source, genes, transcripts, type, start, end, score,
    strand, phase, ID, parent, seqid_str, source_str, attributes;

  seqid = R_splicing_getListElement(pgff, "seqid");
  source = R_splicing_getListElement(pgff, "source");
  genes = R_splicing_getListElement(pgff, "gid");
  transcripts = R_splicing_getListElement(pgff, "tid");
  type = R_splicing_getListElement(pgff, "type");
  start = R_splicing_getListElement(pgff, "start");
  end = R_splicing_getListElement(pgff, "end");
  score = R_splicing_getListElement(pgff, "score");
  strand = R_splicing_getListElement(pgff, "strand");
  phase = R_splicing_getListElement(pgff, "phase");
  ID = R_splicing_getListElement(pgff, "ID");
  parent = R_splicing_getListElement(pgff, "parent");
  seqid_str = R_splicing_getListElement(pgff, "seqid_str");
  source_str = R_splicing_getListElement(pgff, "source_str");
  attributes = R_splicing_getListElement(pgff, "attributes");

  R_splicing_SEXP_to_vector_int_copy(seqid, &gff->seqid);
  R_splicing_SEXP_to_vector_int_copy(source, &gff->source);
  R_splicing_SEXP_to_vector_int_copy(genes, &gff->genes);
  R_splicing_SEXP_to_vector_int_copy(transcripts, &gff->transcripts);
  R_splicing_SEXP_to_vector_int_copy(type, &gff->type);
  R_splicing_SEXP_to_vector_int_copy(start, &gff->start);
  R_splicing_SEXP_to_vector_int_copy(end, &gff->end);
  R_splicing_SEXP_to_vector_copy(score, &gff->score);
  R_splicing_SEXP_to_vector_int_copy(strand, &gff->strand);
  R_splicing_SEXP_to_vector_int_copy(phase, &gff->phase);
  R_splicing_SEXP_to_strvector_copy(ID, &gff->ID);
  R_splicing_SEXP_to_vector_int_copy(parent, &gff->parent);
  R_splicing_SEXP_to_strvector_copy(seqid_str, &gff->seqids);
  R_splicing_SEXP_to_strvector_copy(source_str, &gff->sources);

  gff->n = GET_LENGTH(start);
  gff->nogenes = GET_LENGTH(genes);
  gff->notranscripts = GET_LENGTH(transcripts);

  return 0;
}

int R_splicing_SEXP_to_exonset(SEXP pexons, splicing_exonset_t *exons) {
  
  SEXP seqid_str, seqid, start, end;
  
  seqid = R_splicing_getListElement(pexons, "seqid");
  start = R_splicing_getListElement(pexons, "start");
  end = R_splicing_getListElement(pexons, "end");
  seqid_str = R_splicing_getListElement(pexons, "seqid_str");

  R_splicing_SEXP_to_vector_int(seqid, &exons->seqid);
  R_splicing_SEXP_to_vector_int(start, &exons->start);
  R_splicing_SEXP_to_vector_int(end, &exons->end);
  R_splicing_SEXP_to_strvector(seqid_str, &exons->seqids);

  return 0;
}

int R_splicing_SEXP_to_exons(SEXP pexons, splicing_vector_int_t *exons) {
  int i, p, noexons=GET_LENGTH(pexons);  
  int *tmp=(int*) R_alloc(noexons * 2, sizeof(int));
  splicing_vector_int_view(exons, tmp, noexons*2);

  for (i=0, p=0; i<noexons; i++) {
    int *v=INTEGER(VECTOR_ELT(pexons, i));
    VECTOR(*exons)[p++] = v[0];
    VECTOR(*exons)[p++] = v[1];
  }
  
  return 0;
}

int R_splicing_SEXP_to_isoforms(SEXP piso, splicing_vector_int_t *iso) {
  int i, p, veclen=0, noiso=GET_LENGTH(piso);
  int *tmp;
  
  for (i=0; i<noiso; i++) {
    veclen += GET_LENGTH(VECTOR_ELT(piso, i));
  }
  
  tmp = (int*) R_alloc(veclen+noiso, sizeof(int));
  splicing_vector_int_view(iso, tmp, veclen+noiso);
  
  for (i=0, p=0; i<noiso; i++) {
    int j, ilen=GET_LENGTH(VECTOR_ELT(piso, i));
    int *v=INTEGER(VECTOR_ELT(piso, i));
    for (j=0; j<ilen; j++) {
      VECTOR(*iso)[p++] = v[j]-1;
    }
    VECTOR(*iso)[p++] = -1;
  }
  
  return 0;
}

SEXP R_splicing_strvector_to_SEXP(const splicing_strvector_t *strv) {
  int i, n=splicing_strvector_size(strv);
  SEXP res;
  
  PROTECT(res=NEW_CHARACTER(n));
  for (i=0; i<n; i++) {
    SET_STRING_ELT(res, i, 
		   mkChar(splicing_strvector_get(strv, i)));
  }
  
  UNPROTECT(1);
  return res;
}

SEXP R_splicing_gff_to_SEXP(splicing_gff_t *gff) {
  SEXP result, class, names;

  PROTECT(result=NEW_LIST(15));
  SET_VECTOR_ELT(result, 0, R_splicing_strvector_to_SEXP(&gff->seqids));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&gff->seqid));
  SET_VECTOR_ELT(result, 2, R_splicing_strvector_to_SEXP(&gff->sources));
  SET_VECTOR_ELT(result, 3, R_splicing_vector_int_to_SEXP(&gff->source));
  SET_VECTOR_ELT(result, 4, R_splicing_vector_int_to_SEXP(&gff->type));
  SET_VECTOR_ELT(result, 5, R_splicing_vector_int_to_SEXP(&gff->start));
  SET_VECTOR_ELT(result, 6, R_splicing_vector_int_to_SEXP(&gff->end));
  SET_VECTOR_ELT(result, 7, R_splicing_vector_to_SEXP(&gff->score));
  SET_VECTOR_ELT(result, 8, R_splicing_vector_int_to_SEXP(&gff->strand));
  SET_VECTOR_ELT(result, 9, R_splicing_vector_int_to_SEXP(&gff->phase));
  SET_VECTOR_ELT(result, 10, R_NilValue); /* TODO: keep attributes */
  SET_VECTOR_ELT(result, 11, R_splicing_vector_int_to_SEXP(&gff->genes));
  SET_VECTOR_ELT(result, 12,
		 R_splicing_vector_int_to_SEXP(&gff->transcripts));
  SET_VECTOR_ELT(result, 13, R_splicing_strvector_to_SEXP(&gff->ID));
  SET_VECTOR_ELT(result, 14, R_splicing_vector_int_to_SEXP(&gff->parent));

  PROTECT(names=NEW_CHARACTER(15));
  SET_STRING_ELT(names, 0, mkChar("seqid_str"));
  SET_STRING_ELT(names, 1, mkChar("seqid"));
  SET_STRING_ELT(names, 2, mkChar("source_str"));  
  SET_STRING_ELT(names, 3, mkChar("source"));
  SET_STRING_ELT(names, 4, mkChar("type"));
  SET_STRING_ELT(names, 5, mkChar("start"));
  SET_STRING_ELT(names, 6, mkChar("end"));
  SET_STRING_ELT(names, 7, mkChar("score"));
  SET_STRING_ELT(names, 8, mkChar("strand"));
  SET_STRING_ELT(names, 9, mkChar("phase"));
  SET_STRING_ELT(names, 10, mkChar("attributes"));
  SET_STRING_ELT(names, 11, mkChar("gid"));
  SET_STRING_ELT(names, 12, mkChar("tid"));
  SET_STRING_ELT(names, 13, mkChar("ID"));
  SET_STRING_ELT(names, 14, mkChar("parent"));
  SET_NAMES(result, names);

  PROTECT(class=NEW_CHARACTER(1));
  SET_STRING_ELT(class, 0, mkChar("gff3"));
  SET_CLASS(result, class);

  UNPROTECT(3);
  return result;
}

SEXP R_splicing_exonset_to_SEXP(splicing_exonset_t *exons) {
  SEXP result, class, names;
  
  PROTECT(result=NEW_LIST(4));
  SET_VECTOR_ELT(result, 0, R_splicing_strvector_to_SEXP(&exons->seqids));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&exons->seqid));
  SET_VECTOR_ELT(result, 2, R_splicing_vector_int_to_SEXP(&exons->start));
  SET_VECTOR_ELT(result, 3, R_splicing_vector_int_to_SEXP(&exons->end));
  
  PROTECT(names=NEW_CHARACTER(4));
  SET_STRING_ELT(names, 0, mkChar("seqid_str"));
  SET_STRING_ELT(names, 1, mkChar("seqid"));
  SET_STRING_ELT(names, 2, mkChar("start"));
  SET_STRING_ELT(names, 3, mkChar("end"));
  SET_NAMES(result, names);
  
  PROTECT(class=NEW_CHARACTER(1));
  SET_STRING_ELT(class, 0, mkChar("splicingExonset"));
  SET_CLASS(result, class);
  
  UNPROTECT(3);
  return result;
}

SEXP R_splicing_delbit(SEXP pv, SEXP pbit) {
  int *v=INTEGER(pv);
  int bit=INTEGER(pbit)[0]-1;
  int mask=~(1L << bit);
  int i, n=GET_LENGTH(pv);
  SEXP result;
  int *res;

  R_splicing_begin();
  
  PROTECT(result=NEW_INTEGER(n));
  res=INTEGER(result);
  
  for (i=0; i<n; i++) {
    res[i]=v[i] & mask;
  }

  R_splicing_end();

  UNPROTECT(1);
  return result;
}

SEXP R_splicing_getbit(SEXP pv, SEXP pbit) {
  int *v=INTEGER(pv);
  int bit=INTEGER(pbit)[0]-1;
  int i, n=GET_LENGTH(pv);
  SEXP result;
  int *res;
  
  R_splicing_begin();
  
  PROTECT(result=NEW_INTEGER(n));
  res=INTEGER(result);
  
  for (i=0; i<n; i++) {
    res[i]=(v[i] >> bit) & 1;
  }
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
} 

SEXP R_splicing_read_gff(SEXP file) {
  const char *filename=CHAR(STRING_ELT(file, 0));
  FILE *input=fopen(filename, "r");
  splicing_gff_t gff;
  SEXP result;

  R_splicing_begin();
  
  splicing_gff_init(&gff, 50);
  splicing_gff_read(input, &gff);
  PROTECT(result = R_splicing_gff_to_SEXP(&gff));
  splicing_gff_destroy(&gff);
  fclose(input);
  
  R_splicing_end();

  UNPROTECT(1);
  return result;
}

SEXP R_splicing_write_gff(SEXP pgff, SEXP file) {
  const char *filename=CHAR(STRING_ELT(file, 0));
  FILE *output=fopen(filename, "w");
  splicing_gff_t gff;
  SEXP result;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_gff_write(output, &gff);
  fclose(output);

  R_splicing_end();
  
  return R_NilValue;
}

SEXP R_splicing_matchIso(SEXP pgff, SEXP pgene, SEXP pposition, 
			 SEXP pcigar, SEXP poverhang, SEXP preadlength) {

  int i, noreads=GET_LENGTH(pcigar);
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t exstart, exend, position;
  splicing_gff_t gff;
  int overhang=INTEGER(poverhang)[0];
  int readlength=INTEGER(preadlength)[0];
  splicing_matrix_t res;
  SEXP result;
  const char **cigarstr;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  splicing_matrix_init(&res, 0, 0);
  cigarstr = (const char**) R_alloc(noreads, sizeof(char*));
  
  for (i=0; i<noreads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(pcigar, i));
  }

  splicing_matchIso(&gff, gene, &position, cigarstr, overhang, readlength, 
		    &res);

  PROTECT(result = R_splicing_matrix_to_SEXP(&res));

  splicing_matrix_destroy(&res);

  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_parse_cigar(SEXP cigar, SEXP pmaxreadlength) {
  size_t i, n=GET_LENGTH(cigar);
  int maxreadlength=INTEGER(pmaxreadlength)[0];
  SEXP result, names; 
  splicing_vector_int_t numcigar, cigaridx, cigarlength;
  const char **strcigar;
  
  R_splicing_begin();
  
  PROTECT(result = NEW_LIST(3));
  PROTECT(names = NEW_CHARACTER(3));
  SET_STRING_ELT(names, 0, mkChar("cigar"));
  SET_STRING_ELT(names, 1, mkChar("cigaridx"));
  SET_STRING_ELT(names, 2, mkChar("cigarlength"));
  SET_NAMES(result, names);
  
  splicing_vector_int_init(&numcigar, 0);
  splicing_vector_int_init(&cigaridx, 0);
  splicing_vector_int_init(&cigarlength, 0);
  strcigar = (const char**) R_alloc(n, sizeof(char*));

  for (i=0; i<n; i++) {
    strcigar[i] = CHAR(STRING_ELT(cigar, i));
  }
  
  splicing_parse_cigar(strcigar, n, &numcigar, &cigaridx, &cigarlength,
		       maxreadlength);
  
  SET_VECTOR_ELT(result, 0, R_splicing_vector_int_to_SEXP(&numcigar));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&cigaridx));
  SET_VECTOR_ELT(result, 2, R_splicing_vector_int_to_SEXP(&cigarlength));
  
  splicing_vector_int_destroy(&numcigar);
  splicing_vector_int_destroy(&cigaridx);
  splicing_vector_int_destroy(&cigarlength);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_order_matches(SEXP pmatches) {
  SEXP result; 
  splicing_vector_int_t order;
  splicing_matrix_t matches;

  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix(pmatches, &matches);
  splicing_vector_int_init(&order, 0); 
  
  splicing_order_matches(&matches, &order);
  
  PROTECT(result=R_splicing_vector_int_to_SEXP(&order));
  
  splicing_vector_int_destroy(&order);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_miso(SEXP pgff, SEXP pgene, SEXP preads, SEXP preadLength, 
		     SEXP pnochains, SEXP pnoIterations, 
		     SEXP pmaxIterations, SEXP pnoBurnIn, 
		     SEXP pnoLag, SEXP phyperp, SEXP poverhang, SEXP pstart,
		     SEXP palgo, SEXP pstart_psi, SEXP pstop) {
  
  size_t gene=INTEGER(pgene)[0]-1;
  SEXP result, names, class;
  splicing_matrix_t samples;
  splicing_vector_t logLik;
  splicing_gff_t gff;
  splicing_vector_int_t position;
  const char **cigarstr;
  SEXP cigar=R_splicing_getListElement(preads, "cigar");
  int readLength=INTEGER(preadLength)[0];
  int noChains=INTEGER(pnochains)[0];
  int noIterations=INTEGER(pnoIterations)[0];
  int maxIterations=INTEGER(pmaxIterations)[0];
  int noBurnIn=INTEGER(pnoBurnIn)[0];
  int noLag=INTEGER(pnoLag)[0];
  splicing_vector_t hyperp;
  int overhang=INTEGER(poverhang)[0];
  int i, noReads=GET_LENGTH(cigar);
  splicing_matrix_t match_matrix;
  splicing_matrix_t class_templates;
  splicing_vector_t class_counts;
  splicing_miso_rundata_t rundata;
  int start=INTEGER(pstart)[0];
  int algo=INTEGER(palgo)[0];
  splicing_matrix_t start_psi;
  int stop=INTEGER(pstop)[0];

  R_splicing_begin();
  
  splicing_matrix_init(&samples, 0, 0);
  splicing_vector_init(&logLik, 0);
  splicing_matrix_init(&match_matrix, 0, 0);
  splicing_matrix_init(&class_templates, 0, 0);
  splicing_vector_init(&class_counts, 0);

  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(R_splicing_getListElement(preads, "position"),
				&position);
  cigarstr = (const char**) R_alloc(noReads, sizeof(char*));
  R_splicing_SEXP_to_vector(phyperp, &hyperp);

  for (i=0; i<noReads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(cigar, i));
  }

  if (!isNull(pstart_psi)) {
    R_splicing_SEXP_to_matrix(pstart_psi, &start_psi);
  }

  splicing_miso(&gff, gene, &position, cigarstr, readLength, overhang,
		noChains, noIterations, maxIterations, 
		noBurnIn, noLag, &hyperp, algo, start, stop,
		isNull(pstart_psi) ? 0 : &start_psi,
		&samples, &logLik, &match_matrix, &class_templates, 
		&class_counts, /*assignment=*/ 0, &rundata);
  
  PROTECT(result=NEW_LIST(6));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&samples));
  splicing_matrix_destroy(&samples);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&logLik));
  splicing_vector_destroy(&logLik);
  SET_VECTOR_ELT(result, 2, R_splicing_matrix_to_SEXP(&match_matrix));
  splicing_matrix_destroy(&match_matrix);
  SET_VECTOR_ELT(result, 3, R_splicing_matrix_to_SEXP(&class_templates));
  splicing_matrix_destroy(&class_templates);
  SET_VECTOR_ELT(result, 4, R_splicing_vector_to_SEXP(&class_counts));
  splicing_vector_destroy(&class_counts);
  SET_VECTOR_ELT(result, 5, R_splicing_miso_rundata_to_SEXP(&rundata));

  PROTECT(names=NEW_CHARACTER(6));
  SET_STRING_ELT(names, 0, mkChar("samples"));
  SET_STRING_ELT(names, 1, mkChar("logLik"));
  SET_STRING_ELT(names, 2, mkChar("matchMatrix"));
  SET_STRING_ELT(names, 3, mkChar("classTemplates"));
  SET_STRING_ELT(names, 4, mkChar("classCounts"));
  SET_STRING_ELT(names, 5, mkChar("runData"));
  SET_NAMES(result, names);

  PROTECT(class=ScalarString(mkChar("MISOresult")));
  SET_CLASS(result, class);

  R_splicing_end();
  
  UNPROTECT(3);
  return result;
}

SEXP R_splicing_miso_paired(SEXP pgff, SEXP pgene, SEXP preads, 
			    SEXP preadLength, SEXP pnochains, 
			    SEXP pnoIterations, SEXP pmaxIterations,
			    SEXP pnoBurnIn, SEXP pnoLag, SEXP phyperp, 
			    SEXP pfragmentProb, SEXP pfragmentStart, 
			    SEXP pnormalMean, SEXP pnormalVar, 
			    SEXP pnumDevs, SEXP poverhang, SEXP pstart,
			    SEXP pstartpsi,
			    SEXP pstopcond) {
  
  size_t gene=INTEGER(pgene)[0]-1;
  SEXP result, names, class;
  splicing_matrix_t samples;
  splicing_vector_t logLik;
  splicing_gff_t gff;
  splicing_vector_int_t position;
  const char **cigarstr;
  SEXP cigar=R_splicing_getListElement(preads, "cigar");
  int readLength=INTEGER(preadLength)[0];
  int nochains=INTEGER(pnochains)[0];
  int overHang=INTEGER(poverhang)[0];
  int noIterations=INTEGER(pnoIterations)[0];
  int maxIterations=INTEGER(pmaxIterations)[0];
  int noBurnIn=INTEGER(pnoBurnIn)[0];
  int noLag=INTEGER(pnoLag)[0];
  splicing_vector_t hyperp;
  splicing_vector_t fragmentProb;
  int fragmentStart=INTEGER(pfragmentStart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  int i, noReads=GET_LENGTH(cigar);
  splicing_matrix_t match_matrix;
  splicing_matrix_t class_templates;
  splicing_vector_t class_counts;
  splicing_miso_rundata_t rundata;
  int start=INTEGER(pstart)[0];
  splicing_matrix_t startPsi;
  int stopCond=INTEGER(pstopcond)[0];

  R_splicing_begin();
  
  splicing_matrix_init(&samples, 0, 0);
  splicing_vector_init(&logLik, 0);
  splicing_matrix_init(&match_matrix, 0, 0);
  splicing_matrix_init(&class_templates, 0, 0);
  splicing_vector_init(&class_counts, 0);

  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(R_splicing_getListElement(preads, "position"),
				&position);
  cigarstr = (const char**) R_alloc(noReads, sizeof(char*));
  R_splicing_SEXP_to_vector(phyperp, &hyperp);

  for (i=0; i<noReads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(cigar, i));
  }
  
  if (!isNull(pfragmentProb)) {
    R_splicing_SEXP_to_vector(pfragmentProb, &fragmentProb);
  }

  if (!isNull(pstartpsi)) {
    R_splicing_SEXP_to_matrix(pstartpsi, &startPsi);
  }

  splicing_miso_paired(&gff, gene, &position, cigarstr, readLength,
		       overHang, nochains, noIterations, 
		       maxIterations, noBurnIn, noLag, 
		       &hyperp, start, stopCond, 
		       isNull(pstartpsi) ? 0 : &startPsi, 
		       isNull(pfragmentProb) ? 0 : &fragmentProb, 
		       fragmentStart, normalMean, normalVar, numDevs,
		       &samples, &logLik, &match_matrix, &class_templates,
		       &class_counts, /*bin_class_templates=*/ 0, 
		       /*bin_class_counts=*/ 0, /*assignment=*/ 0,
		       &rundata);
  
  PROTECT(result=NEW_LIST(6));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&samples));
  splicing_matrix_destroy(&samples);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&logLik));
  splicing_vector_destroy(&logLik);
  SET_VECTOR_ELT(result, 2, R_splicing_matrix_to_SEXP(&match_matrix));
  splicing_matrix_destroy(&match_matrix);
  SET_VECTOR_ELT(result, 3, R_splicing_matrix_to_SEXP(&class_templates));
  splicing_matrix_destroy(&class_templates);
  SET_VECTOR_ELT(result, 4, R_splicing_vector_to_SEXP(&class_counts));
  splicing_vector_destroy(&class_counts);
  SET_VECTOR_ELT(result, 5, R_splicing_miso_rundata_to_SEXP(&rundata));

  PROTECT(names=NEW_CHARACTER(6));
  SET_STRING_ELT(names, 0, mkChar("samples"));
  SET_STRING_ELT(names, 1, mkChar("logLik"));
  SET_STRING_ELT(names, 2, mkChar("matchMatrix"));
  SET_STRING_ELT(names, 3, mkChar("classTemplates"));
  SET_STRING_ELT(names, 4, mkChar("classCounts"));
  SET_STRING_ELT(names, 5, mkChar("runData"));
  SET_NAMES(result, names);

  PROTECT(class=ScalarString(mkChar("MISOresult")));
  SET_CLASS(result, class);

  R_splicing_end();
  
  UNPROTECT(3);
  return result;
}

SEXP R_splicing_reassign_samples(SEXP pmatches, SEXP pmatch_order, 
				 SEXP ppsi, SEXP pnoiso, SEXP pnochains) {
  SEXP result; 
  splicing_matrix_t matches;
  splicing_vector_int_t match_order;
  splicing_matrix_t psi;
  int noiso=INTEGER(pnoiso)[0];
  int nochains=INTEGER(pnochains)[0];
  splicing_matrix_int_t cresult;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix(pmatches, &matches);
  R_splicing_SEXP_to_vector_int(pmatch_order, &match_order);
  R_splicing_SEXP_to_matrix(ppsi, &psi);
  splicing_matrix_int_init(&cresult, 0, 0);

  splicing_reassign_samples(&matches, &match_order, &psi, noiso, nochains,
			    &cresult);

  PROTECT(result=R_splicing_matrix_int_to_SEXP(&cresult));
  
  splicing_matrix_int_destroy(&cresult);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_reassign_samples_paired(SEXP pmatches, SEXP pmatch_order,
					SEXP ppsi, SEXP pnoiso, 
					SEXP pnochains,
					SEXP pfragmentStart) {

  SEXP result; 
  splicing_matrix_t matches;
  splicing_vector_int_t match_order;
  splicing_matrix_t psi;
  int noiso=INTEGER(pnoiso)[0];
  int nochains=INTEGER(pnochains)[0];
  int fragmentstart=INTEGER(pfragmentStart)[0];
  splicing_matrix_int_t cresult;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix(pmatches, &matches);
  R_splicing_SEXP_to_vector_int(pmatch_order, &match_order);
  R_splicing_SEXP_to_matrix(ppsi, &psi);
  splicing_matrix_int_init(&cresult, 0, 0);

  splicing_reassign_samples_paired(&matches, &match_order, &psi, noiso, 
				   nochains, fragmentstart, &cresult);

  PROTECT(result=R_splicing_matrix_int_to_SEXP(&cresult));
  
  splicing_matrix_int_destroy(&cresult);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}  
  
SEXP R_splicing_mvplogisnorm(SEXP ptheta, SEXP pmu, SEXP psigma, 
			     SEXP plen) {
  SEXP result;
  splicing_vector_t theta, mu;
  double sigma=REAL(psigma)[0];
  int len=INTEGER(plen)[0];
  double score;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(ptheta, &theta);
  R_splicing_SEXP_to_vector(pmu, &mu);
  
  splicing_mvplogisnorm(&theta, &mu, sigma, len, &score);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0] = sigma;
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_score_iso(SEXP ppsi, SEXP pnoiso, SEXP passignment, 
			  SEXP pnoreads, SEXP peffisolen) {
  SEXP result;
  splicing_vector_t psi;
  int noiso=INTEGER(pnoiso)[0];
  splicing_vector_int_t assignment;
  int noreads=INTEGER(pnoreads)[0];
  splicing_vector_int_t effisolen;
  double res;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(ppsi, &psi);
  R_splicing_SEXP_to_vector_int(passignment, &assignment);
  R_splicing_SEXP_to_vector_int(peffisolen, &effisolen);
  
  splicing_score_iso(&psi, noiso, &assignment, noreads, &effisolen, &res);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0] = res;
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_score_iso_paired(SEXP ppsi, SEXP pnoiso, SEXP passignment,
				 SEXP ppisolen, SEXP passscores) {
  SEXP result;
  splicing_vector_t psi;
  int noiso=INTEGER(pnoiso)[0];
  splicing_vector_int_t assignment;
  splicing_vector_int_t pisolen;
  splicing_vector_t assscores;
  double res;

  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(ppsi, &psi);
  R_splicing_SEXP_to_vector_int(passignment, &assignment);
  R_splicing_SEXP_to_vector_int(ppisolen, &pisolen);
  R_splicing_SEXP_to_vector(passscores, &assscores);
  
  splicing_score_iso_paired(&psi, noiso, &assignment, &pisolen, &assscores,
			    &res);

  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0] = res;
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_ldirichlet(SEXP px, SEXP palpha, SEXP plen) {
  SEXP result;
  splicing_vector_t x, alpha;
  int len=INTEGER(plen)[0];
  double res;

  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(px, &x);
  R_splicing_SEXP_to_vector(palpha, &alpha);
 
  splicing_ldirichlet(&x, &alpha, len, &res);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0] = res;
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_mvrnorm(SEXP pmu, SEXP psigma, SEXP plen) {
  SEXP result;
  splicing_matrix_t mu, resalpha;
  double sigma=REAL(psigma)[0];
  int len=INTEGER(plen)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix(pmu, &mu);
  splicing_matrix_init(&resalpha, 0, 0);

  splicing_mvrnorm(&mu, sigma, &resalpha, len);

  PROTECT(result=R_splicing_matrix_to_SEXP(&resalpha));
  
  splicing_matrix_destroy(&resalpha);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_logit_inv(SEXP px, SEXP plen, SEXP pnochains) {
  SEXP result;
  splicing_matrix_t x;
  splicing_matrix_t res;
  int len=INTEGER(plen)[0];
  int nochains=INTEGER(pnochains)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix(px, &x);
  splicing_matrix_init(&res, len, nochains);
  
  splicing_logit_inv(&x, &res, len, nochains);
  
  PROTECT(result=R_splicing_matrix_to_SEXP(&res));
  
  splicing_matrix_destroy(&res);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_score_joint(SEXP passignment, SEXP pnoreads, SEXP pnochains,
			    SEXP ppsi, SEXP phyper, SEXP peffisolen,
			    SEXP pisoscores, SEXP palgorithm,
			    SEXP pmatchmatrix, SEXP passignmentmatrix,
			    SEXP pmatchvector) {
  
  SEXP result;
  splicing_algorithm_t algorithm=INTEGER(palgorithm)[0];
  splicing_matrix_int_t assignment;
  splicing_vector_int_t effisolen;
  splicing_matrix_t psi;
  splicing_vector_t hyper, isoscores;
  splicing_matrix_t matchMatrix, assignmentMatrix;
  splicing_vector_t matchVector;
  int noreads=INTEGER(pnoreads)[0];
  int nochains=INTEGER(pnochains)[0];
  splicing_vector_t score;
  
  R_splicing_begin();

  if (!isNull(passignment)) {
    R_splicing_SEXP_to_matrix_int(passignment, &assignment);
  }
  R_splicing_SEXP_to_vector_int(peffisolen, &effisolen);
  R_splicing_SEXP_to_matrix(ppsi, &psi);
  R_splicing_SEXP_to_vector(phyper, &hyper);
  R_splicing_SEXP_to_vector(pisoscores, &isoscores);

  if (!isNull(pmatchmatrix)) {
    R_splicing_SEXP_to_matrix(pmatchmatrix, &matchMatrix);
  }
  if (!isNull(passignmentmatrix)) {
    R_splicing_SEXP_to_matrix(passignmentmatrix, &assignmentMatrix);
  }
  if (!isNull(pmatchvector)) {
    R_splicing_SEXP_to_vector(pmatchvector, &matchVector);
  }

  splicing_vector_init(&score, 0);

  splicing_score_joint(algorithm, isNull(passignment) ? 0 : &assignment,
		       noreads, nochains, &psi, &hyper,
		       &effisolen, &isoscores,
		       isNull(pmatchmatrix) ? 0 : &matchMatrix,
		       isNull(passignmentmatrix) ? 0 : &assignmentMatrix,
		       isNull(pmatchvector) ? 0 : &matchVector, &score);
  
  PROTECT(result=R_splicing_vector_to_SEXP(&score));
  splicing_vector_destroy(&score);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_score_joint_paired(SEXP passignment, 
				   SEXP pnoreads, SEXP pnochains, SEXP ppsi, 
				   SEXP phyper, SEXP pisolen, 
				   SEXP pisoscores, SEXP passscores, 
				   SEXP pfragmentLength, 
				   SEXP pfragmentStart) {
  
  SEXP result;
  splicing_matrix_int_t assignment;
  int noreads=INTEGER(pnoreads)[0];
  int nochains=INTEGER(pnochains)[0];
  splicing_matrix_t psi;
  splicing_vector_t hyper;
  splicing_vector_int_t isolen;
  splicing_matrix_t isoscores;
  splicing_vector_t assscores;
  splicing_matrix_int_t fragmentLength;
  int fragmentStart=INTEGER(pfragmentStart)[0];
  splicing_vector_t score;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix_int(passignment, &assignment);
  R_splicing_SEXP_to_matrix(ppsi, &psi);
  R_splicing_SEXP_to_vector(phyper, &hyper);
  R_splicing_SEXP_to_vector_int(pisolen, &isolen);
  R_splicing_SEXP_to_matrix(pisoscores, &isoscores);
  R_splicing_SEXP_to_vector(passscores, &assscores);
  R_splicing_SEXP_to_matrix_int(pfragmentLength, &fragmentLength);

  splicing_vector_init(&score, 0);

  splicing_score_joint_paired(&assignment, noreads, nochains, &psi, 
			      &hyper, &isolen, &isoscores, &assscores, 
			      &fragmentLength, fragmentStart, &score);
  
  PROTECT(result=R_splicing_vector_to_SEXP(&score));
  splicing_vector_destroy(&score);

  R_splicing_end();  

  UNPROTECT(1);
  return result;
}

SEXP R_splicing_drift_proposal(SEXP pmode, SEXP ppsi, SEXP palpha, 
			       SEXP psigma, SEXP potherpsi, SEXP potheralpha,
			       SEXP pnoiso, SEXP pnochains, SEXP pstart,
			       SEXP pstart_psi) {
  int mode=INTEGER(pmode)[0];
  splicing_matrix_t psi, alpha, otherpsi, otheralpha, respsi, resalpha;
  double sigma=REAL(psigma)[0];
  double ressigma;
  splicing_vector_t resscore;
  int noiso=INTEGER(pnoiso)[0];
  int nochains=INTEGER(pnochains)[0];
  SEXP result, names;
  int start=INTEGER(pstart)[0];
  splicing_matrix_t start_psi;

  R_splicing_begin();
  
  switch (mode) {

  case 0:
    splicing_matrix_init(&respsi, 0, 0);
    splicing_matrix_init(&resalpha, 0, 0);
    if (!isNull(pstart_psi)) {
      R_splicing_SEXP_to_matrix(pstart_psi, &start_psi);
    }
    splicing_drift_proposal_init(noiso, nochains, &respsi,
			    &resalpha, &ressigma, start, 
			    isNull(pstart_psi) ? 0 : &start_psi,
			    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    PROTECT(result=NEW_LIST(3));
    SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&respsi));
    SET_VECTOR_ELT(result, 1, R_splicing_matrix_to_SEXP(&resalpha));
    SET_VECTOR_ELT(result, 2, ScalarReal(ressigma));
    PROTECT(names=NEW_CHARACTER(3));
    SET_STRING_ELT(names, 0, mkChar("psi"));
    SET_STRING_ELT(names, 1, mkChar("alpha"));
    SET_STRING_ELT(names, 2, mkChar("sigma"));
    SET_NAMES(result, names);
    splicing_matrix_destroy(&respsi);
    splicing_matrix_destroy(&resalpha);
    UNPROTECT(2);
    break;

  case 1:
    R_splicing_SEXP_to_matrix(ppsi, &psi);
    R_splicing_SEXP_to_matrix(palpha, &alpha);
    splicing_matrix_init(&respsi, 0, 0);
    splicing_matrix_init(&resalpha, 0, 0);
    splicing_drift_proposal_propose(noiso, nochains, &alpha, sigma, 
				    &respsi, &resalpha);
    PROTECT(result=NEW_LIST(2));
    SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&respsi));
    SET_VECTOR_ELT(result, 1, R_splicing_matrix_to_SEXP(&resalpha));
    PROTECT(names=NEW_CHARACTER(2));
    SET_STRING_ELT(names, 0, mkChar("psi"));
    SET_STRING_ELT(names, 1, mkChar("alpha"));
    SET_NAMES(result, names);
    splicing_matrix_destroy(&respsi);
    splicing_matrix_destroy(&resalpha);
    UNPROTECT(2);
    break;

  case 2:
    R_splicing_SEXP_to_matrix(ppsi, &psi);
    R_splicing_SEXP_to_matrix(palpha, &alpha);
    R_splicing_SEXP_to_matrix(potherpsi, &otherpsi);
    R_splicing_SEXP_to_matrix(potheralpha, &otheralpha);
    splicing_vector_init(&resscore, 0);
    splicing_drift_proposal_score(noiso, nochains, &psi, &otheralpha,
				  sigma, &resscore);
    PROTECT(result=R_splicing_vector_to_SEXP(&resscore));
    splicing_vector_destroy(&resscore);
    UNPROTECT(1);
    break;
  }

  R_splicing_end();
  
  return result;
}

SEXP R_splicing_metropolis_hastings_ratio(SEXP passignment, SEXP pnoreads,
					  SEXP pnochains, 
					  SEXP ppsiNew, SEXP palphaNew,
					  SEXP ppsi, SEXP palpha, 
					  SEXP psigma, SEXP pnoiso, 
					  SEXP peffisolen, SEXP phyperp,
					  SEXP pisoscores, SEXP pfull) {

  SEXP result, names; 
  splicing_matrix_int_t assignment;
  splicing_vector_int_t effisolen;
  splicing_matrix_t psiNew, alphaNew, psi, alpha;
  splicing_vector_t hyperp, isoscores;
  int noreads=INTEGER(pnoreads)[0]; 
  int nochains=INTEGER(pnochains)[0];
  int noiso=INTEGER(pnoiso)[0];
  double sigma=REAL(psigma)[0];
  int full=INTEGER(pfull)[0];
  splicing_vector_t acceptP, pcJS, ppJS;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_matrix_int(passignment, &assignment);
  R_splicing_SEXP_to_matrix(ppsiNew, &psiNew);
  R_splicing_SEXP_to_matrix(palphaNew, &alphaNew);
  R_splicing_SEXP_to_matrix(ppsi, &psi);
  R_splicing_SEXP_to_matrix(palpha, &alpha);
  R_splicing_SEXP_to_vector_int(peffisolen, &effisolen);
  R_splicing_SEXP_to_vector(phyperp, &hyperp);
  R_splicing_SEXP_to_vector(pisoscores, &isoscores);

  splicing_vector_init(&acceptP, 0);
  splicing_vector_init(&pcJS, 0);
  splicing_vector_init(&ppJS, 0);

  splicing_metropolis_hastings_ratio(SPLICING_ALGO_REASSIGN,
				     &assignment, noreads, nochains, &psiNew, 
				     &alphaNew, &psi, &alpha, sigma, noiso,
				     /*match=*/ 0, /*assignment=*/ 0,
				     /*matches=*/ 0, &effisolen, &hyperp,
				     &isoscores, full, &acceptP, &pcJS, &ppJS);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_to_SEXP(&acceptP));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&pcJS));
  SET_VECTOR_ELT(result, 2, R_splicing_vector_to_SEXP(&ppJS));
  PROTECT(names=NEW_CHARACTER(3));
  SET_STRING_ELT(names, 0, mkChar("acceptP"));
  SET_STRING_ELT(names, 1, mkChar("pcJS"));
  SET_STRING_ELT(names, 2, mkChar("ppJS"));
  SET_NAMES(result, names);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;  
}

SEXP R_splicing_metropolis_hastings_ratio_paired(SEXP pass, 
						 SEXP pnoreads, 
						 SEXP pnochains, 
						 SEXP ppsiNew,
						 SEXP palphaNew, SEXP ppsi, 
						 SEXP palpha, SEXP psigma,
						 SEXP pnoiso, SEXP pisolen,
						 SEXP phyperp, 
						 SEXP pisoscores,
						 SEXP passscores, 
						 SEXP pfragmentLength,
						 SEXP pfragmentStart, 
						 SEXP pfull) {
  SEXP result, names;
  splicing_matrix_int_t ass;
  int noreads=INTEGER(pnoreads)[0];
  int nochains=INTEGER(pnochains)[0];
  splicing_matrix_t psiNew, alphaNew, psi, alpha;
  splicing_vector_t hyperp, assscores;
  double sigma=REAL(psigma)[0];
  int noiso=INTEGER(pnoiso)[0];
  splicing_vector_int_t isolen;
  splicing_matrix_t isoscores;
  splicing_matrix_int_t fragmentLength;
  int fragmentStart=INTEGER(pfragmentStart)[0];
  int full=INTEGER(pfull)[0];
  splicing_vector_t acceptP, pcJS, ppJS;
  
  R_splicing_begin();

  R_splicing_SEXP_to_matrix_int(pass, &ass);
  R_splicing_SEXP_to_matrix(ppsiNew, &psiNew);
  R_splicing_SEXP_to_matrix(palphaNew, &alphaNew);
  R_splicing_SEXP_to_matrix(ppsi, &psi);
  R_splicing_SEXP_to_matrix(palpha, &alpha);
  R_splicing_SEXP_to_vector_int(pisolen, &isolen);
  R_splicing_SEXP_to_vector(phyperp, &hyperp);
  R_splicing_SEXP_to_matrix(pisoscores, &isoscores);
  R_splicing_SEXP_to_vector(passscores, &assscores);
  R_splicing_SEXP_to_matrix_int(pfragmentLength, &fragmentLength);

  splicing_vector_init(&acceptP, 0);
  splicing_vector_init(&pcJS, 0);
  splicing_vector_init(&ppJS, 0);
  
  splicing_metropolis_hastings_ratio_paired(&ass, noreads, nochains, 
					    &psiNew, &alphaNew, &psi,
					    &alpha, sigma, noiso, &isolen,
					    &hyperp, &isoscores, &assscores,
					    &fragmentLength, fragmentStart,
					    full, &acceptP, &pcJS, &ppJS);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_to_SEXP(&acceptP));
  splicing_vector_destroy(&acceptP);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&pcJS));
  splicing_vector_destroy(&pcJS);
  SET_VECTOR_ELT(result, 2, R_splicing_vector_to_SEXP(&ppJS));
  splicing_vector_destroy(&ppJS);
  PROTECT(names=NEW_CHARACTER(3));
  SET_STRING_ELT(names, 0, mkChar("acceptP"));
  SET_STRING_ELT(names, 1, mkChar("pcJS"));
  SET_STRING_ELT(names, 2, mkChar("ppJS"));
  SET_NAMES(result, names);

  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_gff_exon_start_end(SEXP pgff, SEXP pgene) {
  SEXP result, names;
  int gene=INTEGER(pgene)[0]-1;
  splicing_gff_t gff;
  splicing_vector_int_t start, end, idx;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&start, 0);
  splicing_vector_int_init(&end, 0);
  splicing_vector_int_init(&idx, 0);
  
  splicing_gff_exon_start_end(&gff, &start, &end, &idx, gene);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_int_to_SEXP(&start));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&end));
  SET_VECTOR_ELT(result, 2, R_splicing_vector_int_to_SEXP(&idx));
  PROTECT(names=NEW_CHARACTER(3));
  SET_STRING_ELT(names, 0, mkChar("start"));
  SET_STRING_ELT(names, 1, mkChar("end"));
  SET_STRING_ELT(names, 2, mkChar("idx"));
  SET_NAMES(result, names);

  splicing_vector_int_destroy(&idx);
  splicing_vector_int_destroy(&end);
  splicing_vector_int_destroy(&start);

  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}
  
					  
SEXP R_splicing_numeric_cigar(SEXP pexons, SEXP pnoiso, SEXP pgenestart) {

  SEXP result;
  SEXP pexstart=R_splicing_getListElement(pexons, "start");
  SEXP pexend=R_splicing_getListElement(pexons, "end");
  SEXP pexidx=R_splicing_getListElement(pexons, "idx");
  splicing_vector_int_t exstart, exend, exidx, res;
  int noiso=INTEGER(pnoiso)[0];
  int genestart=INTEGER(pgenestart)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector_int(pexstart, &exstart);
  R_splicing_SEXP_to_vector_int(pexend, &exend);  
  R_splicing_SEXP_to_vector_int(pexidx, &exidx);  
  splicing_vector_int_init(&res, 0);

  splicing_numeric_cigar(&exstart, &exend, &exidx, noiso, genestart, &res);
  
  PROTECT(result=R_splicing_vector_int_to_SEXP(&res));
  splicing_vector_int_destroy(&res);

  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_assignment_matrix(SEXP pgff, SEXP pgene, SEXP preadlength, 
				  SEXP poverhang) {
  SEXP result;
  splicing_matrix_t matrix;
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readlength=INTEGER(preadlength)[0];
  int overhang=INTEGER(poverhang)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_matrix_init(&matrix, 0, 0);
  
  splicing_assignment_matrix(&gff, gene, readlength, overhang, &matrix);
  
  PROTECT(result=R_splicing_matrix_to_SEXP(&matrix));
  
  splicing_matrix_destroy(&matrix);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_assignment_matrix_bias(SEXP pgff, SEXP pgene, 
				       SEXP preadlength, SEXP poverhang, 
				       SEXP pbias) {
  SEXP result;
  splicing_matrix_t matrix;
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readlength=INTEGER(preadlength)[0];
  int overhang=INTEGER(poverhang)[0];
  int bias=INTEGER(pbias)[0];
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_matrix_init(&matrix, 0, 0);
  
  splicing_assignment_matrix_bias(&gff, gene, readlength, overhang, bias, 
				  &matrix);
  
  PROTECT(result=R_splicing_matrix_to_SEXP(&matrix));
  
  splicing_matrix_destroy(&matrix);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_create_gene(SEXP pexons, SEXP pisoforms, SEXP pid, 
			    SEXP pseqid, SEXP psource, SEXP pstrand) {
  SEXP result;
  splicing_vector_int_t exons, isoforms;
  const char *id=CHAR(STRING_ELT(pid, 0));
  const char *seqid=CHAR(STRING_ELT(pseqid, 0));
  const char *source=CHAR(STRING_ELT(psource, 0));
  int strand=INTEGER(pstrand)[0];
  splicing_gff_t gff;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_exons(pexons, &exons);
  R_splicing_SEXP_to_isoforms(pisoforms, &isoforms);
  splicing_gff_init(&gff, 0);
  
  splicing_create_gene(&exons, &isoforms, id, seqid, source, strand, &gff);
  
  PROTECT(result=R_splicing_gff_to_SEXP(&gff));
  
  splicing_gff_destroy(&gff);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_simulate_reads(SEXP pgff, SEXP pgene, SEXP pexpression,
			       SEXP pnoreads, SEXP preadLength) {
  SEXP result, names;
  size_t gene=INTEGER(pgene)[0]-1;
  int noreads=INTEGER(pnoreads)[0];
  int readlength=INTEGER(preadLength)[0];
  splicing_gff_t gff;
  splicing_vector_t expression;
  splicing_vector_int_t isoform, position;
  splicing_strvector_t cigar;
  splicing_vector_t sample_prob;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector(pexpression, &expression);
  splicing_vector_int_init(&isoform, 0);
  splicing_vector_int_init(&position, 0);
  splicing_strvector_init(&cigar, 0);
  splicing_vector_init(&sample_prob, 0);

  splicing_simulate_reads(&gff, gene, &expression, noreads, readlength,
			  &isoform, &position, &cigar, &sample_prob);

  PROTECT(result=NEW_LIST(4));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_int_to_SEXP(&isoform));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&position));
  SET_VECTOR_ELT(result, 2, R_splicing_strvector_to_SEXP(&cigar));
  SET_VECTOR_ELT(result, 3, R_splicing_vector_to_SEXP(&sample_prob));
  PROTECT(names=NEW_CHARACTER(4));
  SET_STRING_ELT(names, 0, mkChar("isoform"));
  SET_STRING_ELT(names, 1, mkChar("position"));
  SET_STRING_ELT(names, 2, mkChar("cigar"));
  SET_STRING_ELT(names, 3, mkChar("sampleProb"));
  SET_NAMES(result, names);
  
  splicing_vector_destroy(&sample_prob);
  splicing_strvector_destroy(&cigar);
  splicing_vector_int_destroy(&position);
  splicing_vector_int_destroy(&isoform);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_simulate_paired_reads(SEXP pgff, SEXP pgene, SEXP pexpression,
				      SEXP pnoReads, SEXP preadLength,
				      SEXP pfragmentProb, SEXP pfragmentStart, 
				      SEXP pnormalMean, SEXP pnormalVar,
				      SEXP pnumDevs) {
  
  SEXP result, names;
  size_t gene=INTEGER(pgene)[0]-1;
  int noreads=INTEGER(pnoReads)[0];
  int readLength=INTEGER(preadLength)[0];
  splicing_gff_t gff;
  splicing_vector_t expression, fragmentProb;
  int fragmentStart=INTEGER(pfragmentStart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  splicing_vector_int_t isoform, position;
  splicing_strvector_t cigar;
  splicing_vector_t sampleprob;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector(pexpression, &expression);
  if (!isNull(pfragmentProb)) { 
    R_splicing_SEXP_to_vector(pfragmentProb, &fragmentProb);
  }
  splicing_vector_int_init(&isoform, 0);
  splicing_vector_int_init(&position, 0);
  splicing_strvector_init(&cigar, 0);
  splicing_vector_init(&sampleprob, 0);
  
  splicing_simulate_paired_reads(&gff, gene, &expression, noreads, 
				 readLength, 
				 isNull(pfragmentProb) ? 0 : &fragmentProb, 
				 fragmentStart, normalMean, normalVar, 
				 numDevs, &isoform, &position, &cigar, 
				 &sampleprob);
  
  PROTECT(result=NEW_LIST(4));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_int_to_SEXP(&isoform));
  SET_VECTOR_ELT(result, 1, R_splicing_vector_int_to_SEXP(&position));
  SET_VECTOR_ELT(result, 2, R_splicing_strvector_to_SEXP(&cigar));
  SET_VECTOR_ELT(result, 3, R_splicing_vector_to_SEXP(&sampleprob));
  PROTECT(names=NEW_CHARACTER(4));
  SET_STRING_ELT(names, 0, mkChar("isoform"));
  SET_STRING_ELT(names, 1, mkChar("position"));
  SET_STRING_ELT(names, 2, mkChar("cigar"));
  SET_STRING_ELT(names, 3, mkChar("sampleProb"));
  SET_NAMES(result, names);
  
  splicing_strvector_destroy(&cigar);
  splicing_vector_int_destroy(&position);
  splicing_vector_int_destroy(&isoform);
  splicing_vector_destroy(&sampleprob);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_gff_noiso(SEXP pgff) {
  SEXP result;
  splicing_gff_t gff;
  splicing_vector_int_t noiso;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&noiso, 0);
  
  splicing_gff_noiso(&gff, &noiso);
  
  PROTECT(result=R_splicing_vector_int_to_SEXP(&noiso));
  
  splicing_vector_int_destroy(&noiso);
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_gff_isolength(SEXP pgff) {
  SEXP result;
  splicing_gff_t gff;
  splicing_vector_int_t isolength, isolength_idx;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&isolength, 0);
  splicing_vector_int_init(&isolength_idx, 0);
  
  splicing_gff_isolength(&gff, &isolength, &isolength_idx);
  
  PROTECT(result=R_splicing_vector_int_index_to_SEXP(&isolength, 
						     &isolength_idx));
  
  splicing_vector_int_destroy(&isolength);
  splicing_vector_int_destroy(&isolength_idx);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_solve_gene(SEXP pgff, SEXP pgene, SEXP preadLength, 
			   SEXP pposition, SEXP pcigar, SEXP poverhang, 
			   SEXP pscale) {

  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readLength=INTEGER(preadLength)[0];
  splicing_vector_int_t position;
  const char **cigarstr;
  int overhang=INTEGER(poverhang)[0];
  int scale=LOGICAL(pscale)[0];
  int i, noReads=GET_LENGTH(pcigar);
  splicing_matrix_t match_matrix, assignment_matrix;
  splicing_vector_t nomatch, expression, residuals;
  double rnorm;
  SEXP result, names;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  splicing_matrix_init(&match_matrix, 0, 0);
  splicing_matrix_init(&assignment_matrix, 0, 0);
  splicing_vector_init(&nomatch, 0);
  splicing_vector_init(&expression, 0);
  splicing_vector_init(&residuals, 0);
  cigarstr = (const char**) R_alloc(noReads, sizeof(char*));
  for (i=0; i<noReads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(pcigar, i));
  }
  
  splicing_solve_gene(&gff, gene, readLength, overhang, &position, cigarstr, 
		      &match_matrix, &nomatch, &assignment_matrix, 
		      &expression, &residuals, scale);
  
  PROTECT(result=NEW_LIST(5));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&match_matrix));
  splicing_matrix_destroy(&match_matrix);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&nomatch));
  splicing_vector_destroy(&nomatch);
  SET_VECTOR_ELT(result, 2, R_splicing_matrix_to_SEXP(&assignment_matrix));
  splicing_matrix_destroy(&assignment_matrix);
  SET_VECTOR_ELT(result, 3, R_splicing_vector_to_SEXP(&expression));
  splicing_vector_destroy(&expression);
  SET_VECTOR_ELT(result, 4, R_splicing_vector_to_SEXP(&residuals));
  splicing_vector_destroy(&residuals);
  
  PROTECT(names=NEW_CHARACTER(5));
  SET_STRING_ELT(names, 0, mkChar("match"));
  SET_STRING_ELT(names, 1, mkChar("nomatch"));
  SET_STRING_ELT(names, 2, mkChar("assignment"));
  SET_STRING_ELT(names, 3, mkChar("expression"));
  SET_STRING_ELT(names, 4, mkChar("residuals"));
  SET_NAMES(result, names); 
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;		 
}

SEXP R_splicing_solve_gene_paired(SEXP pgff, SEXP pgene, SEXP preadLength, 
				  SEXP poverhang, SEXP pscale, 
				  SEXP pposition, SEXP pcigar, SEXP pfast,
				  SEXP pfragmentprob, SEXP pfragmentstart, 
				  SEXP pnormalMean, SEXP pnormalVar, 
				  SEXP pnumDevs) {

  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readLength=INTEGER(preadLength)[0];
  int overhang=INTEGER(poverhang)[0];
  int scale=LOGICAL(pscale)[0];
  splicing_vector_int_t position;
  const char **cigarstr;
  int fast=LOGICAL(pfast)[0];
  int i, noReads=GET_LENGTH(pcigar);
  splicing_vector_t fragmentprob;
  int fragmentstart=INTEGER(pfragmentstart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  splicing_matrix_t match_matrix, assignment_matrix;
  splicing_vector_t expression, residuals, nomatch;
  double rnorm;
  SEXP result, names;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  if (!isNull(pfragmentprob)) { 
    R_splicing_SEXP_to_vector(pfragmentprob, &fragmentprob);
  }
  splicing_matrix_init(&match_matrix, 0, 0);
  splicing_matrix_init(&assignment_matrix, 0, 0);
  splicing_vector_init(&nomatch, 0);
  splicing_vector_init(&expression, 0);
  splicing_vector_init(&residuals, 0);
  cigarstr = (const char**) R_alloc(noReads, sizeof(char*));
  for (i=0; i<noReads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(pcigar, i));
  }
  
  splicing_solve_gene_paired(&gff, gene, readLength, overhang, &position,
			     cigarstr, fast,
			     isNull(pfragmentprob) ? 0 : &fragmentprob, 
			     fragmentstart, normalMean, normalVar, numDevs,
			     &match_matrix, &nomatch, &assignment_matrix, 
			     &expression, &residuals, scale);
  
  PROTECT(result=NEW_LIST(5));
  SET_VECTOR_ELT(result, 0, R_splicing_matrix_to_SEXP(&match_matrix));
  splicing_matrix_destroy(&match_matrix);
  SET_VECTOR_ELT(result, 1, R_splicing_vector_to_SEXP(&nomatch));
  splicing_vector_destroy(&nomatch);
  SET_VECTOR_ELT(result, 2, R_splicing_matrix_to_SEXP(&assignment_matrix));
  splicing_matrix_destroy(&assignment_matrix);
  SET_VECTOR_ELT(result, 3, R_splicing_vector_to_SEXP(&expression));
  splicing_vector_destroy(&expression);
  SET_VECTOR_ELT(result, 4, R_splicing_vector_to_SEXP(&residuals));
  splicing_vector_destroy(&residuals);
  
  PROTECT(names=NEW_CHARACTER(5));
  SET_STRING_ELT(names, 0, mkChar("match"));
  SET_STRING_ELT(names, 1, mkChar("nomatch"));
  SET_STRING_ELT(names, 2, mkChar("assignment"));
  SET_STRING_ELT(names, 3, mkChar("expression"));
  SET_STRING_ELT(names, 4, mkChar("residuals"));
  SET_NAMES(result, names); 
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;		 
}

SEXP R_splicing_iso_to_genomic(SEXP pgff, SEXP pgene, SEXP pisoform,
			       SEXP pposition) {
  
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t isoform, position, newposition;
  SEXP result;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pisoform, &isoform);

  PROTECT(result=NEW_INTEGER(GET_LENGTH(pposition)));
  memcpy(INTEGER(result), INTEGER(pposition), 
	 sizeof(int) * GET_LENGTH(pposition));
  splicing_vector_int_view(&newposition, INTEGER(result), 
			   GET_LENGTH(pposition));
  
  splicing_iso_to_genomic(&gff, gene, &isoform, 0, &newposition);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_genomic_to_iso(SEXP pgff, SEXP pgene, SEXP pposition) {

  SEXP result;
  
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t position;
  splicing_matrix_int_t isopos;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  splicing_matrix_int_init(&isopos, 0, 0);

  splicing_genomic_to_iso(&gff, gene, &position, 0, &isopos);
  
  PROTECT(result=R_splicing_matrix_int_to_SEXP(&isopos));
  
  splicing_matrix_int_destroy(&isopos);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_matchIso_paired(SEXP pgff, SEXP pgene, SEXP pposition,
				SEXP pcigar, SEXP preadlength, 
				SEXP poverhang,
				SEXP pfragmentprob, SEXP pfragmentstart, 
				SEXP pnormalMean, SEXP pnormalVar,
				SEXP pnumDevs) {
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t position;
  size_t i, noreads=GET_LENGTH(pposition);
  const char **cigarstr;
  int readlength=INTEGER(preadlength)[0];
  int overhang=INTEGER(poverhang)[0];
  splicing_vector_t fragmentprob;
  int fragmentstart=INTEGER(pfragmentstart)[0];
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  splicing_matrix_t result;
  splicing_matrix_int_t fragmentLength;
  SEXP rresult;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  R_splicing_SEXP_to_vector_int(pposition, &position);
  cigarstr = (const char**) R_alloc(noreads, sizeof(char*));
  for (i=0; i<noreads; i++) {
    cigarstr[i] = CHAR(STRING_ELT(pcigar, i));
  }
  if (!isNull(pfragmentprob)) {
    R_splicing_SEXP_to_vector(pfragmentprob, &fragmentprob);
  }
  splicing_matrix_init(&result, 0, 0);
  splicing_matrix_int_init(&fragmentLength, 0, 0);
  
  splicing_matchIso_paired(&gff, gene, &position, cigarstr, readlength, 
			   overhang, 
			   isNull(pfragmentprob) ? 0 : &fragmentprob, 
			   fragmentstart, normalMean, normalVar, numDevs,
			   &result, &fragmentLength);

  PROTECT(rresult=NEW_LIST(2));
  SET_VECTOR_ELT(rresult, 0, R_splicing_matrix_to_SEXP(&result));
  SET_VECTOR_ELT(rresult, 1, R_splicing_matrix_int_to_SEXP(&fragmentLength));

  splicing_matrix_destroy(&result);
  splicing_matrix_int_destroy(&fragmentLength);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return rresult;
}

SEXP R_splicing_paired_assignment_matrix(SEXP pgff, SEXP pgene, 
					 SEXP preadlength, SEXP poverhang,
					 SEXP pbias, SEXP pfast,
					 SEXP pfragmentprob,
					 SEXP pfragmentstart, 
					 SEXP pnormalMean, SEXP pnormalVar,
					 SEXP pnumDevs) {
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readlength=INTEGER(preadlength)[0];
  int overhang=INTEGER(poverhang)[0];
  int bias=INTEGER(pbias)[0];
  splicing_vector_t fragmentprob;
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];  
  int fragmentstart=INTEGER(pfragmentstart)[0];
  int fast=LOGICAL(pfast)[0];
  splicing_matrix_t matrix;
  SEXP result;

  R_splicing_begin();
  
  R_splicing_SEXP_to_gff(pgff, &gff);
  if (!isNull(pfragmentprob)) {
    R_splicing_SEXP_to_vector(pfragmentprob, &fragmentprob);
  }
  splicing_matrix_init(&matrix, 0, 0);

  if (!fast) {
    splicing_paired_assignment_matrix_bias(&gff, gene, readlength, overhang, 
			   bias, isNull(pfragmentprob) ? 0 : &fragmentprob,
			   fragmentstart, normalMean, normalVar,
			   numDevs, &matrix);
  } else {
    splicing_paired_assignment_matrix_fast_bias(&gff, gene, readlength, 
	       overhang, bias, isNull(pfragmentprob) ? 0 : &fragmentprob,
	       fragmentstart, normalMean, normalVar, numDevs, &matrix);
  }
  
  PROTECT(result=R_splicing_matrix_to_SEXP(&matrix));
  
  splicing_matrix_destroy(&matrix);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_normal_fragment(SEXP pnormalMean, SEXP pnormalVar, 
				SEXP pnumDevs, SEXP pminLength) { 
  SEXP result, names;
  double normalMean=REAL(pnormalMean)[0];
  double normalVar=REAL(pnormalVar)[0];
  double numDevs=REAL(pnumDevs)[0];
  int minLength=INTEGER(pminLength)[0];
  splicing_vector_t fragmentProb;
  int fragmentStart;
  
  R_splicing_begin();
  
  splicing_vector_init(&fragmentProb, 0);
  splicing_normal_fragment(normalMean, normalVar, numDevs, minLength, 
			   &fragmentProb, &fragmentStart);

  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, R_splicing_vector_to_SEXP(&fragmentProb));
  splicing_vector_destroy(&fragmentProb);
  SET_VECTOR_ELT(result, 1, ScalarInteger(fragmentStart));
  PROTECT(names=NEW_CHARACTER(2));
  SET_STRING_ELT(names, 0, mkChar("fragmentProb"));
  SET_STRING_ELT(names, 1, mkChar("fragmentStart"));
  SET_NAMES(result, names);
  
  R_splicing_end();
  
  UNPROTECT(2);
  return result;
}

SEXP R_splicing_gene_complexity(SEXP pgff, SEXP pgene, 
				SEXP preadlength, 
				SEXP poverhang, SEXP ptype, 
				SEXP pnorm, SEXP ppaired, SEXP pfast, 
				SEXP pfragmentprob, SEXP pfragmentstart, 
				SEXP pnormalmean, SEXP pnormalvar, 
				SEXP pnumdevs) {
  SEXP result;
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int readlength=INTEGER(preadlength)[0];
  int overhang=INTEGER(poverhang)[0];
  splicing_complexity_t type=INTEGER(ptype)[0];
  splicing_norm_t norm=INTEGER(pnorm)[0];
  int paired=LOGICAL(ppaired)[0];
  int fast=LOGICAL(pfast)[0];
  splicing_vector_t fragmentprob;
  int fragmentstart=INTEGER(pfragmentstart)[0];
  double normalmean=REAL(pnormalmean)[0];
  double normalvar=REAL(pnormalvar)[0];
  double numdevs=REAL(pnumdevs)[0];
  double complexity;

  R_splicing_begin();

  R_splicing_SEXP_to_gff(pgff, &gff);
  if (!isNull(pfragmentprob)) {
    R_splicing_SEXP_to_vector(pfragmentprob, &fragmentprob);
  }
  
  splicing_gene_complexity(&gff, gene, readlength, overhang, type, norm, 
			   paired, fast, 
			   isNull(pfragmentprob) ? 0 : &fragmentprob,
			   fragmentstart, normalmean, normalvar, numdevs,
			   &complexity);
  
  PROTECT(result=ScalarReal(complexity));

  R_splicing_end();

  UNPROTECT(1);  
  return result;
}

SEXP R_splicing_read_sambam(SEXP pfilename) {
  const char *filename=CHAR(STRING_ELT(pfilename, 0));
  splicing_reads_t reads;
  SEXP result;
  
  R_splicing_begin();

  splicing_reads_init(&reads);
  splicing_read_sambam(filename, SPLICING_SAMBAM_AUTO, &reads);
  PROTECT(result=R_splicing_reads_to_SEXP(&reads));
  splicing_reads_destroy(&reads);
  
  R_splicing_end();

  UNPROTECT(1);
  return result;
}

SEXP R_splicing_read_sambam_region(SEXP pfilename, SEXP pregion) {
  const char *filename=CHAR(STRING_ELT(pfilename, 0));
  const char *region=CHAR(STRING_ELT(pregion, 0));
  splicing_reads_t reads;
  SEXP result;
  
  R_splicing_begin();
  
  splicing_reads_init(&reads);
  splicing_read_sambam_region(filename, /*indexfile=*/ 0, 
			      SPLICING_SAMBAM_AUTO, region, &reads);
  PROTECT(result=R_splicing_reads_to_SEXP(&reads));
  splicing_reads_destroy(&reads);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_gtf2gff3(SEXP gtf, SEXP gid, SEXP tid) {
  SEXP result, rownames, colnames;
  SEXP Tseqname, Tsource, Tfeature, Tstart, Tend, Tscore, Tstrand, Tframe;
  SEXP Fseqid, Fsource, Ftype, Fstart, Fend, Fscore, Fstrand, Fphase, 
    Fattributes;
    
  int Tnrow=GET_LENGTH(gid), Fnrow=0;
  int noGenes=1, noTranscripts=1, i=0, j=0;
  const char *acttid, *actgid;
  int gStart, gEnd, tStart, tEnd; 
  double min, max;
  char attr[2000];

  R_splicing_begin();

  /* Input fields */

  Tseqname = R_splicing_getListElement(gtf, "seqname");
  Tsource  = R_splicing_getListElement(gtf, "source");
  Tfeature = R_splicing_getListElement(gtf, "feature");
  Tstart   = R_splicing_getListElement(gtf, "start");
  Tend     = R_splicing_getListElement(gtf, "end");
  Tscore   = R_splicing_getListElement(gtf, "score");
  Tstrand  = R_splicing_getListElement(gtf, "strand");
  Tframe   = R_splicing_getListElement(gtf, "frame");

  /* --------------------------------------------------------------- */
  /* Set up output                                                   */
  /* --------------------------------------------------------------- */

  /* Count number of unique genes, number of unique transcripts
     We use that '(gid, tid)' is sorted, first according to 'gid',
     then according to 'tid' */

  actgid=CHARACTER_VALUE(STRING_ELT(gid, 0));
  acttid=CHARACTER_VALUE(STRING_ELT(tid, 0));
  for (i=1; i<Tnrow; i++) {
    const char *mgid=CHARACTER_VALUE(STRING_ELT(gid, i));
    const char *mtid=CHARACTER_VALUE(STRING_ELT(tid, i));
    if (strcmp(actgid, mgid)) { noGenes++;       actgid=mgid; }
    if (strcmp(acttid, mtid)) { noTranscripts++; acttid=mtid; }
  }
  Fnrow = Tnrow + noGenes + noTranscripts;
  
  /* Allocate memory for the new rows */

  PROTECT(result=NEW_LIST(9));
  SET_VECTOR_ELT(result, 0, NEW_CHARACTER(Fnrow)); /* seqid */
  SET_VECTOR_ELT(result, 1, NEW_CHARACTER(Fnrow)); /* source */
  SET_VECTOR_ELT(result, 2, NEW_CHARACTER(Fnrow)); /* type */
  SET_VECTOR_ELT(result, 3, NEW_INTEGER(Fnrow));   /* start */
  SET_VECTOR_ELT(result, 4, NEW_INTEGER(Fnrow));   /* end */
  SET_VECTOR_ELT(result, 5, NEW_NUMERIC(Fnrow));   /* score */
  SET_VECTOR_ELT(result, 6, NEW_CHARACTER(Fnrow)); /* strand */
  SET_VECTOR_ELT(result, 7, NEW_CHARACTER(Fnrow)); /* phase */
  SET_VECTOR_ELT(result, 8, NEW_CHARACTER(Fnrow)); /* attributes */
  Fseqid      = VECTOR_ELT(result, 0);
  Fsource     = VECTOR_ELT(result, 1);
  Ftype       = VECTOR_ELT(result, 2);
  Fstart      = VECTOR_ELT(result, 3);
  Fend        = VECTOR_ELT(result, 4);
  Fscore      = VECTOR_ELT(result, 5);
  Fstrand     = VECTOR_ELT(result, 6);
  Fphase      = VECTOR_ELT(result, 7);
  Fattributes = VECTOR_ELT(result, 8);

  /* Names */

  PROTECT(colnames=NEW_CHARACTER(9));
  SET_STRING_ELT(colnames, 0, mkChar("seqid"));
  SET_STRING_ELT(colnames, 1, mkChar("source"));
  SET_STRING_ELT(colnames, 2, mkChar("type"));
  SET_STRING_ELT(colnames, 3, mkChar("start"));
  SET_STRING_ELT(colnames, 4, mkChar("end"));
  SET_STRING_ELT(colnames, 5, mkChar("score"));
  SET_STRING_ELT(colnames, 6, mkChar("strand"));
  SET_STRING_ELT(colnames, 7, mkChar("phase"));
  SET_STRING_ELT(colnames, 8, mkChar("attributes"));
  SET_NAMES(result, colnames);

  PROTECT(rownames = NEW_INTEGER(Fnrow));
  for (i=0, j=1; i<Fnrow; i++, j++) { INTEGER(rownames)[i]=j; }
  SET_ATTR(result, install("row.names"), rownames);

  /* Class */

  SET_CLASS(result, ScalarString(CREATE_STRING_VECTOR("data.frame")));  
  
  /* --------------------------------------------------------------- */
  /* Conversion                                                      */
  /* --------------------------------------------------------------- */

  /* We just go over the original data, and add the gene and mRNA
     entries. For this we need to find where each transcript and gene 
     starts and ends. */

#define GID(I) CHARACTER_VALUE(STRING_ELT(gid, (I)))
#define TID(I) CHARACTER_VALUE(STRING_ELT(tid, (I)))

  /* Macro from an abstract read.
     Place the beginning of the gene in 'S', the end in 'E', 
     start searching from 'F', and also place start position in 'min'
     and end position in 'max'. */

#define READ_GENE(S,E,F,MIN,MAX)			                  \
  if ((F) >= Tnrow) {  							  \
    (S)=-1;                                                               \
  } else {                                                                \
    int from=(F);							  \
    (MIN)=INTEGER(Tstart)[from];                                          \
    (MAX)=INTEGER(Tend)[from];				                  \
    (S)=from; (E)=from+1;				                  \
    while ((E) < Tnrow && !strcmp(GID(S), GID(E))) {                      \
      if (INTEGER(Tstart)[(E)] < (MIN)) { (MIN)=INTEGER(Tstart)[(E)]; }   \
      if (INTEGER(Tend)[(E)] > (MAX)) { (MAX)=INTEGER(Tend)[(E)]; }       \
      (E)=(E)+1;					                  \
    }                                                                     \
    (E)=(E)-1;						                  \
  }

#define READ_ISOFORM(S,E,F,O,MIN,MAX)                                     \
  if ((F) > (O)) {                                                        \
    (S)=-1;                                                               \
  } else {                                                                \
    int from=(F);                                                         \
    (S)=from; (E)=from+1;                                                 \
    while ((E) <= (O) && !strcmp(TID(S), TID(E))) {                       \
      (E)=(E)+1;                                                          \
    }                                                                     \
    (E)=(E)-1;                                                            \
  }

  i=0;
  READ_GENE(gStart, gEnd, 0, min, max);
  while (gStart != -1) {

    /* Allow CTRL + C */
    R_CheckUserInterrupt();

    /* Entry for 'gene' */
    SET_STRING_ELT(Fseqid, i, STRING_ELT(Tseqname, gStart));
    SET_STRING_ELT(Fsource, i, STRING_ELT(Tsource, gStart));
    SET_STRING_ELT(Ftype, i, mkChar("gene"));
    INTEGER(Fstart)[i]=min;
    INTEGER(Fend)[i]=max;
    REAL(Fscore)[i] = REAL(Tscore)[gStart];
    SET_STRING_ELT(Fstrand, i, STRING_ELT(Tstrand, gStart));
    SET_STRING_ELT(Fphase, i, STRING_ELT(Tframe, gStart));
    snprintf(attr, sizeof(attr)/sizeof(char)-sizeof(char), 
	     "ID=%s", GID(gStart));
    SET_STRING_ELT(Fattributes, i, mkChar(attr));
    i++;

    /* Do the individual isoforms */
    READ_ISOFORM(tStart, tEnd, gStart, gEnd, min, max);
    while (tStart != -1) {
      int no, nExon=1, nCDS=1, nStartC=1, nStopC=1, nOther=1;
      
      /* Do the isoform */
      SET_STRING_ELT(Fseqid, i, STRING_ELT(Tseqname, tStart));
      SET_STRING_ELT(Fsource, i, STRING_ELT(Tsource, tStart));
      SET_STRING_ELT(Ftype, i, mkChar("mRNA"));
      INTEGER(Fstart)[i]=min;
      INTEGER(Fend)[i]=max;
      REAL(Fscore)[i] = REAL(Tscore)[tStart];
      SET_STRING_ELT(Fstrand, i, STRING_ELT(Tstrand, tStart));
      SET_STRING_ELT(Fphase, i, STRING_ELT(Tframe, tStart));
      snprintf(attr, sizeof(attr)/sizeof(char)-sizeof(char), 
	       "ID=%s;Parent=%s", TID(tStart), GID(tStart));
      SET_STRING_ELT(Fattributes, i, mkChar(attr));
      i++;
      for (j=tStart; j<=tEnd; j++) {
	SET_STRING_ELT(Fseqid, i, STRING_ELT(Tseqname, j));
	SET_STRING_ELT(Fsource, i, STRING_ELT(Tsource, j));
	SET_STRING_ELT(Ftype, i, STRING_ELT(Tfeature, j));
	INTEGER(Fstart)[i]=INTEGER(Tstart)[j];
	INTEGER(Fend)[i]=INTEGER(Tend)[j];
	REAL(Fscore)[i] = REAL(Tscore)[j];
	SET_STRING_ELT(Fstrand, i, STRING_ELT(Tstrand, j));
	SET_STRING_ELT(Fphase, i, STRING_ELT(Tframe, j));
	if (!strcmp(CHARACTER_VALUE(STRING_ELT(Tfeature, j)), "exon")) {
	  no=nExon++;
	} else if (!strcmp(CHARACTER_VALUE(STRING_ELT(Tfeature, j)), 
			   "CDS")) {
	  no=nCDS++;
	} else if (!strcmp(CHARACTER_VALUE(STRING_ELT(Tfeature, j)), 
			   "start_codon")) {
	  no=nStartC++;
	} else if (!strcmp(CHARACTER_VALUE(STRING_ELT(Tfeature, j)), 
			   "stop_codon")) {
	  no=nStopC++;
	} else {
	  no=nOther++;
	}
	snprintf(attr, sizeof(attr)/sizeof(char)-sizeof(char),
		 "ID=%s:%s:%i;Parent=%s", TID(j), 
		 CHARACTER_VALUE(STRING_ELT(Tfeature, j)),
		 no, TID(j));
	SET_STRING_ELT(Fattributes, i, mkChar(attr));
	i++;
      }
      
      READ_ISOFORM(tStart, tEnd, tEnd+1, gEnd, min, max);
    }

    READ_GENE(gStart, gEnd, gEnd+1, min, max);
  }

  for (; i<Fnrow; i++) {
    INTEGER(Fstart)[i] = -1;
  }

  R_splicing_end();

  UNPROTECT(3);
  return result;
}

SEXP R_splicing_noexons_one(SEXP pgff, SEXP pgene) {
  
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  splicing_vector_int_t noexons;
  SEXP result;
  
  R_splicing_begin();

  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&noexons, 0);
  
  splicing_gff_noexons_one(&gff, gene, &noexons);
  
  PROTECT(result=R_splicing_vector_int_to_SEXP(&noexons));
  splicing_vector_int_destroy(&noexons);

  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_writemiso(SEXP pmisoresult, SEXP pfile) {
  SEXP psamples=R_splicing_getListElement(pmisoresult, "samples");
  SEXP plogLik=R_splicing_getListElement(pmisoresult, "logLik");
  SEXP pmatchMatrix=R_splicing_getListElement(pmisoresult, "matchMatrix");
  int noReads=INTEGER(GET_DIM(pmatchMatrix))[1];
  SEXP pclassTemplates=R_splicing_getListElement(pmisoresult, 
						 "classTemplates");
  int noClasses=INTEGER(GET_DIM(pclassTemplates))[1];
  SEXP pclassCounts=R_splicing_getListElement(pmisoresult, "classCounts");
  SEXP prunData=R_splicing_getListElement(pmisoresult, "runData");
  SEXP pgeneStructure=R_splicing_getListElement(pmisoresult, "geneStructure");
  
  double *samples=REAL(psamples);
  int noSamples=INTEGER(GET_DIM(psamples))[1];
  double *logLik=REAL(plogLik);
  double *matchMatrix=REAL(pmatchMatrix);
  double *classTemplates=REAL(pclassTemplates);
  double *classCounts=REAL(pclassCounts);
  int noIso=INTEGER(R_splicing_getListElement(prunData, "noIso"))[0];
  int noIters=INTEGER(R_splicing_getListElement(prunData, "noIters"))[0];
  int noBurnIn=INTEGER(R_splicing_getListElement(prunData, "noBurnIn"))[0];
  int noLag=INTEGER(R_splicing_getListElement(prunData, "noLag"))[0];
  int noAccepted=INTEGER(R_splicing_getListElement(prunData, 
						   "noAccepted"))[0];
  int noRejected=INTEGER(R_splicing_getListElement(prunData, 
						   "noRejected"))[0];
  int noChains=INTEGER(R_splicing_getListElement(prunData, 
						 "noChains"))[0];
  splicing_gff_t gff;
  const char *filename=CHAR(STRING_ELT(pfile, 0));
  FILE *file=fopen(filename, "w");
  int i, j, idx;
  
  R_splicing_begin();

  R_splicing_SEXP_to_gff(pgeneStructure, &gff);
  
  fputs("[runData]\n", file);
  fprintf(file, "noReads: %i\n", noReads);
  fprintf(file, "noIso: %i\n", noIso);
  fprintf(file, "noIters: %i\n", noIters);
  fprintf(file, "noBurnIn: %i\n", noBurnIn);
  fprintf(file, "noLag: %i\n", noLag);
  fprintf(file, "noAccepted: %i\n", noAccepted);
  fprintf(file, "noRejected: %i\n", noRejected);
  fprintf(file, "nochains: %i\n", noChains);
  fprintf(file, "noSamples: %i\n", noSamples);
  
  fputs("\n[geneStructure]\n", file);
  splicing_gff_write(file, &gff);
  
  fputs("\n[classTemplates]\n", file);
  for (i=0, idx=0; i<noClasses; i++) {
    for (j=0; j<noIso-1; j++) {
      fprintf(file, "%g ", classTemplates[idx++]);
    }
    fprintf(file, "%g\n", classTemplates[idx++]);
  }

  fputs("\n[classCounts]\n", file);
  for (i=0, idx=0; i<noClasses; i++) {
    fprintf(file, "%g\n", classCounts[idx++]);
  }
  
  fputs("\n[matchMatrix]\n", file);
  for (i=0, idx=0; i<noReads; i++) { 
    for (j=0; j<noIso-1; j++) {
      fprintf(file, "%g ", matchMatrix[idx++]);
    }
    fprintf(file, "%g\n", matchMatrix[idx++]);
  }

  fputs("\n[samples]\n", file);
  for (i=0, idx=0; i<noSamples; i++) {
    for (j=0; j<noIso-1; j++) {
      fprintf(file, "%g ", samples[idx++]);
    }
    fprintf(file, "%g\n", samples[idx++]);
  }
  
  fputs("\n[logLik]\n", file);
  for (i=0, idx=0; i<noSamples; i++) {
    fprintf(file, "%g\n", logLik[idx++]);
  }

  fclose(file);
  R_splicing_end();  
  return R_NilValue;
}

SEXP R_splicing_writelinear(SEXP plinresult, SEXP pfile) {
  SEXP pmatch=R_splicing_getListElement(plinresult, "match");
  SEXP pnomatch=R_splicing_getListElement(plinresult, "nomatch");
  SEXP passignment=R_splicing_getListElement(plinresult, "assignment");
  SEXP pexpression=R_splicing_getListElement(plinresult, "expression");
  SEXP presiduals=R_splicing_getListElement(plinresult, "residuals");
  SEXP pgeneStructure=R_splicing_getListElement(plinresult, "geneStructure");
  
  int noClasses=INTEGER(GET_DIM(passignment))[1];
  int noReads=INTEGER(GET_DIM(pmatch))[1];
  int noIso=INTEGER(GET_DIM(pmatch))[0];

  double *matchMatrix=REAL(pmatch);
  double *nomatch=REAL(pnomatch);
  double *assignment=REAL(passignment);
  double *expression=REAL(pexpression);
  double *residuals=REAL(presiduals);

  splicing_gff_t gff;

  const char *filename=CHAR(STRING_ELT(pfile, 0));
  FILE *file=fopen(filename, "w");
  int i, j, idx;
  
  R_splicing_begin();

  R_splicing_SEXP_to_gff(pgeneStructure, &gff);

  fputs("[runData]\n", file);
  fprintf(file, "noReads: %i\n", noReads);
  fprintf(file, "noIso: %i\n", noIso);

  fputs("\n[expression]\n", file);
  for (i=0; i<noIso-1; i++) {
    fprintf(file, "%g ", expression[i]);
  } 
  fprintf(file, "%g\n", expression[i]);

  fputs("\n[residuals]\n", file);
  for (i=0; i<noIso-1; i++) {
    fprintf(file, "%g ", residuals[i]);
  } 
  fprintf(file, "%g\n", residuals[i]);

  fputs("\n[geneStructure]\n", file);
  splicing_gff_write(file, &gff);
  
  fputs("\n[classTemplates]\n", file);
  for (i=0, idx=0; i<noClasses; i++) {
    for (j=0; j<noIso-1; j++) {
      fprintf(file, "%i ", assignment[idx++] > 0 ? 1 : 0);
    }
    fprintf(file, "%i\n", assignment[idx++] > 0 ? 1 : 0);
  }

  fputs("\n[classCounts]\n", file);
  for (i=0, idx=0; i<noClasses; i++) {
    fprintf(file, "%g\n", nomatch[idx++]);
  }

  fputs("\n[matchMatrix]\n", file);
  for (i=0, idx=0; i<noReads; i++) { 
    for (j=0; j<noIso-1; j++) {
      fprintf(file, "%g ", matchMatrix[idx++]);
    }
    fprintf(file, "%g\n", matchMatrix[idx++]);
  }  
  
  R_splicing_end();
  return R_NilValue;
}

SEXP R_splicing_sam2bam(SEXP pinfile, SEXP poutfile) {
  const char *infile=CHAR(STRING_ELT(pinfile, 0));
  const char *outfile=CHAR(STRING_ELT(poutfile, 0));

  R_splicing_begin();

  splicing_sam2bam(infile, outfile);

  R_splicing_end();

  return R_NilValue;
}

SEXP R_splicing_bam_sort(SEXP pinfile, SEXP poutprefix, SEXP pkey) {
  
  const char *infile=CHAR(STRING_ELT(pinfile, 0));
  const char *outprefix=CHAR(STRING_ELT(poutprefix, 0));
  int key=INTEGER(pkey)[0];
  
  R_splicing_begin();
  
  splicing_bam_sort(infile, outprefix, key);
  
  R_splicing_end();
  
  return R_NilValue;
}

SEXP R_splicing_bam_index(SEXP pfilename) {
  const char *filename=CHAR(STRING_ELT(pfilename, 0));
  
  R_splicing_begin();
  
  splicing_bam_index(filename);
  
  R_splicing_end();
  
  return R_NilValue;
}

SEXP R_splicing_rng_get_dirichlet(SEXP palpha) {
  splicing_vector_t alpha;
  splicing_vector_t res;
  SEXP result;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_vector(palpha, &alpha);
  splicing_vector_init(&res, 0);
  
  splicing_rng_get_dirichlet(&splicing_rng_default, &alpha, &res);
  
  PROTECT(result=R_splicing_vector_to_SEXP(&res));
  splicing_vector_destroy(&res);

  R_splicing_end();

  UNPROTECT(1);
  return result;
}

SEXP R_splicing_genomic_to_iso_all(SEXP pgff, SEXP pgene, SEXP pposition) {
  SEXP result;
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int position=INTEGER(pposition)[0];
  splicing_vector_int_t cresult;
  
  R_splicing_begin();

  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&cresult, 0);
  
  splicing_genomic_to_iso_all(&gff, gene, position, /*converter=*/ 0, 
			      &cresult);

  PROTECT(result=R_splicing_vector_int_to_SEXP(&cresult));
  splicing_vector_int_destroy(&cresult);

  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_iso_to_genomic_all(SEXP pgff, SEXP pgene, SEXP pposition) {
  SEXP result;
  splicing_gff_t gff;
  int gene=INTEGER(pgene)[0]-1;
  int position=INTEGER(pposition)[0];
  splicing_vector_int_t cresult;
  
  R_splicing_begin();

  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_vector_int_init(&cresult, 0);
  
  splicing_iso_to_genomic_all(&gff, gene, position, /*converter=*/ 0, 
			      &cresult);

  PROTECT(result=R_splicing_vector_int_to_SEXP(&cresult));
  splicing_vector_int_destroy(&cresult);

  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}
				  
SEXP R_splicing_constitutive_exons(SEXP pgff, SEXP pmin_length, SEXP pmode) {
  SEXP result;
  splicing_gff_t gff;
  splicing_exonset_t exons;
  int min_length=INTEGER(pmin_length)[0];
  splicing_constitutive_mode_t mode=INTEGER(pmode)[0];
  
  R_splicing_begin();

  R_splicing_SEXP_to_gff(pgff, &gff);
  splicing_gff_constitutive_exons(&gff, &exons, min_length, mode);
  PROTECT(result = R_splicing_exonset_to_SEXP(&exons));
  splicing_exonset_destroy(&exons);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_estimate_fraglength(SEXP pexons, SEXP preadsfile, 
				    SEXP preads) {
  SEXP result;
  splicing_exonset_t exons;
  const char *readsfile= isNull(preadsfile) ? 0 : 
    CHAR(STRING_ELT(preadsfile, 0));
  splicing_reads_t reads;
  splicing_vector_int_t fraglen;
  
  R_splicing_begin();
  
  R_splicing_SEXP_to_exonset(pexons, &exons);
  if (!isNull(preads)) { R_splicing_SEXP_to_reads(preads, &reads); }
  splicing_vector_int_init(&fraglen, 0);
  splicing_estimate_fragment_length(&exons, readsfile, 
				    isNull(preads) ? 0 : &reads, &fraglen);
  PROTECT(result = R_splicing_vector_int_to_SEXP(&fraglen));
  splicing_vector_int_destroy(&fraglen);
  
  R_splicing_end();
  
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_version() {
  const char *version;
  SEXP result;
  splicing_version(&version, /*major=*/ 0, /*minor=*/ 0, /*patch=*/ 0);
  PROTECT(result=NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(version));
  UNPROTECT(1);
  return result;
}

SEXP R_splicing_reindex_gff(SEXP pgff) {
  splicing_gff_t gff;
  SEXP result;

  R_splicing_begin();

  R_splicing_SEXP_to_gff_copy(pgff, &gff);
  splicing_gff_reindex(&gff);
  PROTECT(result = R_splicing_gff_to_SEXP(&gff));
  splicing_gff_destroy(&gff);

  R_splicing_end();

  UNPROTECT(1);
  return result;
}

SEXP R_splicing_cr_matrix(SEXP p_gff, SEXP p_gene, SEXP p_ass_mat,
			  SEXP p_read_len, SEXP p_overhang, SEXP p_expr,
			  SEXP p_paired, SEXP p_fast_ass_mat,
			  SEXP p_fragment_prob, SEXP p_fragment_start,
			  SEXP p_normal_mean, SEXP p_normal_var,
			  SEXP p_num_devs) {
  SEXP result;
  splicing_gff_t gff;
  size_t gene = INTEGER(p_gene)[0] - 1;
  splicing_matrix_t ass_mat;
  int read_len = INTEGER(p_read_len)[0];
  int overhang = INTEGER(p_overhang)[0];
  splicing_vector_t expr;
  int paired = LOGICAL(p_paired)[0];
  int fast_ass_mat = LOGICAL(p_fast_ass_mat)[0];
  splicing_vector_t fragment_prob;
  int fragment_start = INTEGER(p_fragment_start)[0];
  double normal_mean = REAL(p_normal_mean)[0];
  double normal_var = REAL(p_normal_var)[0];
  double num_devs = REAL(p_num_devs)[0];
  splicing_matrix_t c_result;

  R_splicing_SEXP_to_gff(p_gff, &gff);
  R_splicing_SEXP_to_vector(p_expr, &expr);
  if (!isNull(p_ass_mat)) {
    R_splicing_SEXP_to_matrix(p_ass_mat, &ass_mat);
  }
  if (!isNull(p_fragment_prob)) {
    R_splicing_SEXP_to_vector(p_fragment_prob, &fragment_prob);
  }

  splicing_matrix_init(&c_result, 0, 0);
  splicing_cr_matrix(&gff, gene, isNull(p_ass_mat) ? 0 : &ass_mat,
		     read_len, overhang, &expr, paired, fast_ass_mat,
		     isNull(p_fragment_prob) ? 0 : &fragment_prob,
		     fragment_start, normal_mean, normal_var, num_devs,
		     &c_result);

  PROTECT(result = R_splicing_matrix_to_SEXP(&c_result));
  splicing_matrix_destroy(&c_result);

  UNPROTECT(1);
  return result;
}
