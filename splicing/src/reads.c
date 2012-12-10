
#include <ctype.h>
#include <string.h>
#include <unistd.h>

#include "splicing.h"
#include "splicing_error.h"

#include "sam.h"
#include "Rsplicing.h"

bam_index_t *bam_index_load_core(FILE *fp);

int splicing_reads_init(splicing_reads_t *reads) {
  reads->noPairs = reads->noSingles = 0;
  reads->paired = 0;

  SPLICING_CHECK(splicing_strvector_init(&reads->chrname, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->chrname);
  SPLICING_CHECK(splicing_vector_int_init(&reads->chrlen, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->chrlen);
  SPLICING_CHECK(splicing_vector_int_init(&reads->chr, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->chr);
  SPLICING_CHECK(splicing_strvector_init(&reads->qname, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->qname);
  SPLICING_CHECK(splicing_strvector_init(&reads->cigar, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->cigar);
  SPLICING_CHECK(splicing_vector_int_init(&reads->position, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->position);
  SPLICING_CHECK(splicing_vector_int_init(&reads->flags, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->flags);
  SPLICING_CHECK(splicing_vector_int_init(&reads->pairpos, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->pairpos);
  SPLICING_CHECK(splicing_vector_int_init(&reads->mapq, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->mapq);
  SPLICING_CHECK(splicing_vector_int_init(&reads->rnext, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->rnext);
  SPLICING_CHECK(splicing_vector_int_init(&reads->tlen, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->tlen);
  SPLICING_CHECK(splicing_strvector_init(&reads->seq, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->seq);
  SPLICING_CHECK(splicing_strvector_init(&reads->qual, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->qual);
  SPLICING_CHECK(splicing_vector_int_init(&reads->mypair, 0));
  SPLICING_FINALLY(splicing_vector_int_destroy, &reads->mypair);
  SPLICING_CHECK(splicing_strvector_init(&reads->attributes, 0));
  SPLICING_FINALLY(splicing_strvector_destroy, &reads->attributes);
  
  SPLICING_FINALLY_CLEAN(15);

  return 0;
}

void splicing_reads_destroy(splicing_reads_t *reads) {
  splicing_strvector_destroy(&reads->chrname);
  splicing_vector_int_destroy(&reads->chrlen);
  splicing_vector_int_destroy(&reads->chr);
  splicing_strvector_destroy(&reads->qname);
  splicing_strvector_destroy(&reads->cigar);
  splicing_vector_int_destroy(&reads->position);
  splicing_vector_int_destroy(&reads->flags);
  splicing_vector_int_destroy(&reads->pairpos);
  splicing_vector_int_destroy(&reads->mapq);
  splicing_vector_int_destroy(&reads->rnext);
  splicing_vector_int_destroy(&reads->tlen);
  splicing_strvector_destroy(&reads->seq);
  splicing_strvector_destroy(&reads->qual);
  splicing_vector_int_destroy(&reads->mypair);
  splicing_strvector_destroy(&reads->attributes);
}

int splicing_reads_clear(splicing_reads_t *reads) {
  splicing_strvector_clear(&reads->chrname);
  splicing_vector_int_clear(&reads->chrlen);
  splicing_vector_int_clear(&reads->chr);
  splicing_strvector_clear(&reads->qname);
  splicing_strvector_clear(&reads->cigar);
  splicing_vector_int_clear(&reads->position);
  splicing_vector_int_clear(&reads->flags);
  splicing_vector_int_clear(&reads->pairpos);
  splicing_vector_int_clear(&reads->mapq);
  splicing_vector_int_clear(&reads->rnext);
  splicing_vector_int_clear(&reads->tlen);
  splicing_strvector_clear(&reads->seq);
  splicing_strvector_clear(&reads->qual);
  splicing_vector_int_clear(&reads->mypair);
  splicing_strvector_clear(&reads->attributes);
  reads->noPairs = reads->noSingles = 0;
  return 0;
}

void splicing_i_bam_destroy1(bam1_t *read) {
  bam_destroy1(read);
}

typedef struct splicing_i_read_sambam_sort_data_t {
  splicing_vector_int_t *chr;
  splicing_vector_int_t *pos;
  splicing_vector_int_t *flag;
  splicing_vector_int_t *pairpos;
} splicing_i_read_sambam_sort_data_t;

int splicing_i_cmp_reads(void *pdata, const void *a, const void *b) {
  const splicing_i_read_sambam_sort_data_t *data = 
    (splicing_i_read_sambam_sort_data_t *) pdata;
  int aa = *(int*) a, bb = *(int*) b;
  int apos, bpos, arpos, brpos, aapos, bbpos;

  if (VECTOR(*data->chr)[aa] < VECTOR(*data->chr)[bb]) { 
    return -1; 
  } else if (VECTOR(*data->chr)[aa] > VECTOR(*data->chr)[bb]) { 
    return 1;
  }
  
  aapos=apos=VECTOR(*data->pos)[aa];
  bbpos=bpos=VECTOR(*data->pos)[bb];

  arpos=VECTOR(*data->pairpos)[aa]; 
  if (arpos > apos) { int tmp=arpos; arpos=apos; apos=tmp; } 
  brpos=VECTOR(*data->pairpos)[bb];
  if (brpos > bpos) { int tmp=brpos; brpos=bpos; bpos=tmp; } 

  if (arpos < brpos) { return -1; } else if (arpos > brpos) { return 1; }
  if (apos  < bpos ) { return -1; } else if (apos  > bpos ) { return 1; }
  if (aapos < bbpos) { return -1; } else if (aapos > bbpos) { return 1; }  
  
  return 0;
}

#define REMAINING (sizeof(buffer)-(bufptr-buffer)*sizeof(char))

int splicing_i_add_read(splicing_reads_t *reads, const bam1_t *read) {
  int i, ncigar=read->core.n_cigar;
  char buffer[4096];
  char *cigarcode="MIDNSHP";
  char *bufptr=buffer;
  uint32_t *actcigar=bam1_cigar(read);
  uint8_t *s = bam1_seq(read), *t = bam1_qual(read);
  
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->position,
					       read->core.pos+1));
  SPLICING_CHECK(splicing_strvector_append2(&reads->qname, 
					    bam1_qname(read), 
					    read->core.l_qname-1));
  if (read->core.n_cigar == 0) {
    SPLICING_CHECK(splicing_strvector_append(&reads->cigar, "*"));
  } else {
    for (i=0; i<ncigar; i++) {
      int l=snprintf(bufptr, REMAINING,
		     "%i%c", (int) (actcigar[i] >> BAM_CIGAR_SHIFT), 
		     cigarcode[actcigar[i] & BAM_CIGAR_MASK]);
      bufptr += l;
    }
    SPLICING_CHECK(splicing_strvector_append2(&reads->cigar, buffer,
					      bufptr-buffer));
  }
  
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->chr, read->core.tid));
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->pairpos, 
					       read->core.mpos+1));
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->flags, 
					       read->core.flag));
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->rnext, 
					       read->core.mtid));
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->mapq, 
					       read->core.qual));
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->tlen,
					       read->core.isize));
    
  if (!read->core.l_qseq) {
    SPLICING_CHECK(splicing_strvector_append(&reads->seq, "*"));
  } else {
    bufptr=buffer;
    for (i=0; i<read->core.l_qseq; i++) { 
      *bufptr = bam_nt16_rev_table[bam1_seqi(s, i)];
      bufptr++;
    }
    SPLICING_CHECK(splicing_strvector_append2(&reads->seq, buffer, 
					      bufptr-buffer));
  }
  
  if (!read->core.l_qseq || t[0] == 0xff) {
    SPLICING_CHECK(splicing_strvector_append(&reads->qual, "*"));
  } else {
    bufptr=buffer;
    for (i=0; i<read->core.l_qseq; i++) { 
      *bufptr = t[i] + 33;
      bufptr++;
    }
    SPLICING_CHECK(splicing_strvector_append2(&reads->qual, buffer,
					      bufptr-buffer));
  }

  SPLICING_CHECK(splicing_vector_int_push_back(&reads->mypair, -1));

  /* Extra columns */
  
  bufptr = buffer;
  s = bam1_aux(read);
  i=0;
  while (s < read->data + read->data_len) {
    uint8_t type, key[2];
    int l;
    key[0] = s[0]; key[1] = s[1];
    s += 2; type = *s; ++s;
    if (i!=0) { *bufptr = '\t'; bufptr++; }
    l=snprintf(bufptr, REMAINING, "%c%c:", (int) key[0], (int) key[1]);
    bufptr += l;
    switch (type) { 
      uint8_t sub_type;
      int32_t n;
      int j;
    case 'A': 
      l=snprintf(bufptr, REMAINING, "A:%c", (int) *s); ++s; bufptr += l;
      break;
    case 'C':
      l=snprintf(bufptr, REMAINING, "i:%u", *(uint8_t*)s); ++s; bufptr += l;
      break;
    case 'c':
      l=snprintf(bufptr, REMAINING, "i:%i", *(int8_t*)s); ++s; bufptr += l;
      break;
    case 'S':
      l=snprintf(bufptr, REMAINING, "i:%u", *(uint16_t*)s); s+=2; bufptr += l;
      break;
    case 's': 
      l=snprintf(bufptr, REMAINING, "i:%i", *(int16_t*)s); s+=2; bufptr += l;
      break;
    case 'I':
      l=snprintf(bufptr, REMAINING, "i:%u", *(uint32_t*)s); s+=4; bufptr +=l;
      break;
    case 'i':
      l=snprintf(bufptr, REMAINING, "i:%i", *(int32_t*)s); s+=4; bufptr += l;
      break;
    case 'f':
      l=snprintf(bufptr, REMAINING, "f:%g", *(float*)s); s+=4; bufptr += l;
      break;
    case 'd':
      l=snprintf(bufptr, REMAINING, "d:%lg", *(double*)s); s+=8; bufptr += l;
      break;
    case 'z':
    case 'H':
      l=snprintf(bufptr, REMAINING, "%c:", type); bufptr += l;
      while (*s) { 
	if (REMAINING > 0) { *bufptr = *s; bufptr++; }
	s++;
      }
      break;
    case 'B':
      sub_type = *(s++);
      memcpy(&n, s, 4);
      s += 4;
      l=snprintf(bufptr, REMAINING, "%c:%c", type, sub_type); bufptr += l;
      for (j = 0; j < n; ++j) {
	l=snprintf(bufptr, REMAINING, ","); bufptr += l;
	switch (sub_type) {
	case 'c':
	  l=snprintf(bufptr, REMAINING, "%i", *(int8_t*)s); ++s; bufptr += l;
	  break;
	case 'C': 
	  l=snprintf(bufptr, REMAINING, "%u", *(uint8_t*)s); ++s; bufptr += l;
	  break;
	case 's': 
	  l=snprintf(bufptr, REMAINING, "%i", *(int16_t*)s); 
	  s += 2; bufptr += l;
	  break;
	case 'S':
	  l=snprintf(bufptr, REMAINING, "%u", *(uint16_t*)s); 
	  s += 2; bufptr += l;
	  break;
	case 'i': 
	  l=snprintf(bufptr, REMAINING, "%i", *(int32_t*)s); 
	  s += 4; bufptr += l;
	  break;
	case 'I': 
	  l=snprintf(bufptr, REMAINING, "%u", *(uint32_t*)s); 
	  s += 4; bufptr += l;
	  break;
	case 'f':
	  l=snprintf(bufptr, REMAINING, "%g", *(float*)s); 
	  s += 4; bufptr += l;
	  break;
	}
      }
    }
    
    i++;
  }
  SPLICING_CHECK(splicing_strvector_append2(&reads->attributes, buffer, 
					    bufptr-buffer));
  
  /* Single-end or paired-end? */
  
  if (read->core.mpos < 0) {
    reads->noSingles  += 1;
  } else {
    reads->noPairs += 1;
  }
  
  return 0;
}

int splicing_i_add_read_minimal(splicing_reads_t *reads, const bam1_t *read) {

  int i, ncigar=read->core.n_cigar;
  char buffer[4096];
  char *cigarcode="MIDNSHP";
  char *bufptr=buffer;
  uint32_t *actcigar=bam1_cigar(read);

  /* Reads without a pair are skipped already here */
  if (read->core.mpos < 0) { return 0; }

  SPLICING_CHECK(splicing_vector_int_push_back(&reads->position, 
					       read->core.pos+1));
  SPLICING_CHECK(splicing_strvector_append2(&reads->qname, bam1_qname(read), 
					    read->core.l_qname-1));
  if (read->core.n_cigar == 0) {
    SPLICING_CHECK(splicing_strvector_append(&reads->cigar, "*"));
  } else {
    for (i=0; i<ncigar; i++) {
      int l=snprintf(bufptr, REMAINING,
		     "%i%c", (int) (actcigar[i] >> BAM_CIGAR_SHIFT), 
		     cigarcode[actcigar[i] & BAM_CIGAR_MASK]);
      bufptr += l;
    }
    SPLICING_CHECK(splicing_strvector_append2(&reads->cigar, buffer,
					      bufptr-buffer));
  }
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->chr, read->core.tid));
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->pairpos, 
					       read->core.mpos+1));
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->flags, 
					       read->core.flag));
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->rnext, 
					       read->core.mtid));
  /* MAPQ skipped */
  /* TLEN skipped */
  /* SEQ skipped */
  /* QUAL skipped */
  
  SPLICING_CHECK(splicing_vector_int_push_back(&reads->mypair, -1));
  
  /* ATTRIBUTES skipped */
  
  reads->noPairs += 1;

  return 0;
}

#undef REMAINING

int splicing_i_order_reads(splicing_reads_t *reads) {
  splicing_vector_char_t taken;
  splicing_vector_int_t idx;
  splicing_vector_t idx2;
  size_t i, pos, noReads=splicing_vector_int_size(&reads->position);
  splicing_i_read_sambam_sort_data_t data = 
    { &reads->chr, &reads->position, &reads->flags, &reads->pairpos };
  int nonpaired_mask = ~0x1;

  SPLICING_CHECK(splicing_vector_int_init(&idx, noReads));
  SPLICING_FINALLY(splicing_vector_int_destroy, &idx);
  SPLICING_CHECK(splicing_vector_init(&idx2, noReads));
  SPLICING_FINALLY(splicing_vector_destroy, &idx2);
  SPLICING_CHECK(splicing_vector_char_init(&taken, noReads));
  SPLICING_FINALLY(splicing_vector_char_destroy, &taken);
  for (i=0; i<noReads; i++) { VECTOR(idx)[i] = i; }
  splicing_qsort_r(VECTOR(idx), noReads, sizeof(int), (void*) &data, 
		   splicing_i_cmp_reads);
  
  for (i=0, pos=0; i<noReads; i++) {
    int curr=VECTOR(idx)[i];
    if (VECTOR(taken)[curr]) {
      /* already done, do nothing */
    } else if (VECTOR(reads->pairpos)[curr] == 0) { 
      /* has no pair, add it */
      VECTOR(reads->mypair)[pos] = -1;
      VECTOR(idx2)[pos++] = curr;
    } else if (VECTOR(reads->pairpos)[curr] < VECTOR(reads->position)[curr] 
	       && !VECTOR(taken)[curr])  {
      /* pair is missing */
      VECTOR(reads->mypair)[pos] = -1;
      VECTOR(idx2)[pos++] = curr;
      VECTOR(reads->pairpos)[curr] = 0;
      VECTOR(reads->flags)[curr] &= nonpaired_mask;
      reads->noPairs--;
      reads->noSingles++;
    } else {
      /* search for pair, forward */
      int needle=VECTOR(reads->pairpos)[curr], ppos, found;
      for (ppos=i+1, found=0; ! found && ppos < noReads; ppos++) {
	int rp=VECTOR(idx)[ppos];
	found = VECTOR(reads->position)[rp] == needle && 
	  VECTOR(reads->pairpos)[rp] == VECTOR(reads->position)[curr] &&
	  0x80 & VECTOR(reads->flags)[rp] && ! VECTOR(taken)[rp];
      }
      if (found) { 
	int rp=VECTOR(idx)[ppos-1];
	VECTOR(reads->mypair)[pos] = pos+1;
	VECTOR(idx2)[pos++] = curr;
	VECTOR(reads->mypair)[pos] = pos-1;
	VECTOR(idx2)[pos++] = rp;
	VECTOR(taken)[rp] = 1;
      } else {
	/* pair is missing */
	VECTOR(reads->mypair)[pos] = -1;
	VECTOR(idx2)[pos++] = curr;
	VECTOR(reads->pairpos)[curr] = 0;
	VECTOR(reads->flags)[curr] &= nonpaired_mask;
	reads->noPairs--;
	reads->noSingles++;
      }
    }
  }

  splicing_vector_char_destroy(&taken);
  SPLICING_FINALLY_CLEAN(1);

  /* We have the correct order now, reorder the vectors */
  splicing_vector_int_iindex(&reads->chr, &idx2);
  splicing_strvector_permute(&reads->qname, &idx2);
  splicing_strvector_permute(&reads->cigar, &idx2);
  splicing_vector_int_iindex(&reads->position, &idx2);
  splicing_vector_int_iindex(&reads->flags, &idx2);
  splicing_vector_int_iindex(&reads->pairpos, &idx2);
  splicing_vector_int_iindex(&reads->mapq, &idx2);
  splicing_vector_int_iindex(&reads->rnext, &idx2);
  splicing_vector_int_iindex(&reads->tlen, &idx2);
  splicing_strvector_permute(&reads->seq, &idx2);
  splicing_strvector_permute(&reads->attributes, &idx2);
  
  splicing_vector_destroy(&idx2);
  splicing_vector_int_destroy(&idx);
  SPLICING_FINALLY_CLEAN(2);

  return 0;
}

/* TODO: use qname as well for pairing */

int splicing_read_sambam(const char *filename,
			 splicing_sambam_type_t filetype,
			 splicing_reads_t *reads) {

  samfile_t *infile=0;
  int bytesread;
  bam1_t *read = bam_init1();
  int i;
  char *mode_r="r", *mode_rb="rb", *mode=mode_rb;

  bam_verbose=0;

  SPLICING_FINALLY(splicing_i_bam_destroy1, read);

  splicing_reads_clear(reads);

  switch (filetype) {
    int flen;
  case SPLICING_SAMBAM_AUTO:
    flen=strlen(filename);
    if (flen >= 4 && !strncmp(filename + flen - 4, ".sam", 4)) { 
      mode = mode_r; 
      filetype=SPLICING_SAMBAM_SAM;
    } else {
      filetype=SPLICING_SAMBAM_BAM;
    }
    break;
  case SPLICING_SAMBAM_SAM:
    mode = mode_r;
    break;
  case SPLICING_SAMBAM_BAM:
    mode = mode_rb;
    break;
  }

  infile=samopen(filename, mode, /*aux=*/ 0);
  if (!infile) { 
    SPLICING_ERROR("Cannot open SAM/BAM file", SPLICING_EFILE);
  }

  SPLICING_CHECK(splicing_vector_int_resize(&reads->chrlen,
					    infile->header->n_targets));
  for (i=0; i<infile->header->n_targets; i++) {
    char *s=infile->header->target_name[i];
    SPLICING_CHECK(splicing_strvector_append(&reads->chrname, s));
    VECTOR(reads->chrlen)[i] = infile->header->target_len[i];
  }

  while ( (bytesread=samread(infile, read)) >= 0) {
    SPLICING_CHECK(splicing_i_add_read(reads, read));
  }

  if (bytesread < -1) {
    SPLICING_WARNING("Truncated SAM/BAM file");
  }

  /* We need to find the pairs of the paired-end reads, if there is any.
     For this we first order the reads according to their RNAME
     (chromosome), plus their position. Then we search for the second
     pair of each first pair. */

  if (reads->noPairs != 0) {
    SPLICING_CHECK(splicing_i_order_reads(reads));
  }

  reads->noPairs /= 2;
  reads->paired = reads->noPairs ? 1 : 0;
    
  bam_destroy1(read);
  SPLICING_FINALLY_CLEAN(1);

  samclose(infile);
  
  return 0;
}

int splicing_i_read_sambam_cb(const bam1_t *read, void *data) {
  splicing_reads_t *reads = (splicing_reads_t*) data;
  SPLICING_CHECK(splicing_i_add_read(reads, read));
  return 0;
}

int splicing_read_sambam_region(const char *filename,
				const char *indexfile,
				splicing_sambam_type_t filetype,
				const char *region,
				splicing_reads_t *reads) {

  samfile_t *infile=0;
  bam1_t *read = bam_init1();
  int i;
  char *mode_r="r", *mode_rb="rb", *mode=mode_rb;
  char *myindexfile=(char*) indexfile;
  bam_index_t *idx=0;
  int tid, beg, end, result;

  bam_verbose=0;

  SPLICING_FINALLY(splicing_i_bam_destroy1, read);

  splicing_reads_clear(reads);

  switch (filetype) {
    int flen;
  case SPLICING_SAMBAM_AUTO:
    flen=strlen(filename);
    if (flen >= 4 && !strncmp(filename + flen - 4, ".sam", 4)) { 
      mode = mode_r; 
      filetype=SPLICING_SAMBAM_SAM;
    } else {
      filetype=SPLICING_SAMBAM_BAM;
    }
    break;
  case SPLICING_SAMBAM_SAM:
    mode = mode_r;
    break;
  case SPLICING_SAMBAM_BAM:
    mode = mode_rb;
    break;
  }

  /* Load index, if available */

  if (filetype == SPLICING_SAMBAM_BAM) { 
    FILE * ifp;
    if (!indexfile) { 
      int flen=strlen(filename);
      myindexfile = malloc(flen + 5);
      strcpy(myindexfile, filename);
      strcat(myindexfile, ".bai");
    }
    ifp = fopen(myindexfile, "rb");
    if (!indexfile) { free(myindexfile); }
    if (ifp) { 
      idx = bam_index_load_core(ifp);
      SPLICING_FINALLY(bam_index_destroy, idx);
      fclose(ifp);
    }
    if (!idx) {
      SPLICING_ERROR("Cannot read BAM index file", SPLICING_EFILE);
    }
  }

  infile=samopen(filename, mode, /*aux=*/ 0);
  if (!infile) { 
    SPLICING_ERROR("Cannot open SAM/BAM file", SPLICING_EFILE);
  }

  SPLICING_CHECK(splicing_vector_int_resize(&reads->chrlen,
					    infile->header->n_targets));
  for (i=0; i<infile->header->n_targets; i++) {
    char *s=infile->header->target_name[i];
    SPLICING_CHECK(splicing_strvector_append(&reads->chrname, s));
    VECTOR(reads->chrlen)[i] = infile->header->target_len[i];
  }

  bam_parse_region(infile->header, region, &tid, &beg, &end);
  if (tid < 0) { 
    SPLICING_ERROR("Unknown reference name", SPLICING_EINVAL);
  }
  
  result = bam_fetch(infile->x.bam, idx, tid, beg, end, reads, 
		     splicing_i_read_sambam_cb);
  if (result < 0) {
    SPLICING_ERROR("Truncated BAM file or corrupt BAM index file",
		   SPLICING_EINVAL);
  }

  /* We need to find the pairs of the paired-end reads, if there is any.
     For this we first order the reads according to their RNAME
     (chromosome), plus their position. Then we search for the second
     pair of each first pair. */

  if (reads->noPairs != 0) {
    SPLICING_CHECK(splicing_i_order_reads(reads));
  }

  reads->noPairs /= 2;
  reads->paired = reads->noPairs ? 1 : 0;
  
  if (idx) {
    bam_index_destroy(idx);
    SPLICING_FINALLY_CLEAN(1);
  }
  
  bam_destroy1(read);
  SPLICING_FINALLY_CLEAN(1);

  samclose(infile);  

  return 0;
}

int splicing_sam2bam(const char *infile, const char *outfile) {

  samfile_t *in = 0, *out = 0;
  bam1_t *b = bam_init1();
  int r;

  bam_verbose=0;

  if ((in = samopen(infile, "r", /*aux=*/ 0)) == 0) {
    SPLICING_ERROR("Failed to open sam file for reading.", SPLICING_EFILE);
  }
  if ((out = samopen(outfile, "wbh", in->header)) == 0) {  
    samclose(in);
    SPLICING_ERROR("Failed to open bam file for writing.", SPLICING_EFILE);
  }
  
  while ((r = samread(in, b)) >= 0) {
    samwrite(out, b);
  }
  if (r < -1) {
    SPLICING_WARNING("Truncated file");
  }
  bam_destroy1(b);

  samclose(in);
  samclose(out);
  
  return 0;
}

typedef bam1_t *bam1_p;
extern int g_is_by_qname;
extern void sort_blocks(int n, int k, bam1_p *buf, const char *prefix,
			const bam_header_t *h, int is_stdout);
int bam_merge_core(int by_qname, const char *out, const char *headers, 
		   int n, char * const *fn, int flag, const char *reg);

/* TODO: free memory in case of error */

int splicing_bam_sort(const char *infile, const char *outprefix, 
		      splicing_bam_sort_key_t key) {
  
  int is_by_qname = (key == SPLICING_BAM_KEY_QNAME);
  int is_stdout = 0;
  size_t max_mem = 500000000;
  
  int n, ret, k, i;
  size_t mem;
  bam_header_t *header;
  bamFile fp;
  bam1_t *b, **buf;

  bam_verbose=0;
  
  g_is_by_qname = is_by_qname;
  n = k = 0; mem = 0;
  fp = bam_open(infile, "r");
  if (fp == 0) {
    SPLICING_ERROR("Cannot open BAM file", SPLICING_EFILE);
  }
  header = bam_header_read(fp);
  buf = (bam1_t**)calloc(max_mem / BAM_CORE_SIZE, sizeof(bam1_t*));
  // write sub files
  for (;;) {
    if (buf[k] == 0) buf[k] = (bam1_t*)calloc(1, sizeof(bam1_t));
    b = buf[k];
    if ((ret = bam_read1(fp, b)) < 0) break;
    mem += ret;
    ++k;
    if (mem >= max_mem) {
      sort_blocks(n++, k, buf, outprefix, header, 0);
      mem = 0; k = 0;
    }
  }
  if (ret != -1)
    SPLICING_WARNING("[bam_sort_core] truncated file. Continue anyway.\n");
  if (n == 0) 
    sort_blocks(-1, k, buf, outprefix, header, is_stdout);
  else { // then merge
    char **fns, *fnout;
    fprintf(stderr, "[bam_sort_core] merging from %d files...\n", n+1);
    sort_blocks(n++, k, buf, outprefix, header, 0);
    fnout = (char*)calloc(strlen(outprefix) + 20, 1);
    if (is_stdout) sprintf(fnout, "-");
    else sprintf(fnout, "%s.bam", outprefix);
    fns = (char**)calloc(n, sizeof(char*));
    for (i = 0; i < n; ++i) {
      fns[i] = (char*)calloc(strlen(outprefix) + 20, 1);
      sprintf(fns[i], "%s.%.4d.bam", outprefix, i);
    }
    bam_merge_core(is_by_qname, fnout, 0, n, fns, 0, 0);
    free(fnout);
    for (i = 0; i < n; ++i) {
      unlink(fns[i]);
      free(fns[i]);
    }
    free(fns);
  }
  for (k = 0; k < max_mem / BAM_CORE_SIZE; ++k) {
    if (buf[k]) {
      free(buf[k]->data);
      free(buf[k]);
    }
  }
  free(buf);
  bam_header_destroy(header);
  bam_close(fp);
  
  return 0;
}

int splicing_bam_index(const char *filename) {
  bam_verbose=0;
  bam_index_build(filename);
  return 0;
}

int splicing_i_estfraglen_cb(const bam1_t *read, void *data) {
  splicing_reads_t *reads = (splicing_reads_t*) data;
  SPLICING_CHECK(splicing_i_add_read_minimal(reads, read));
  return 0;  
}

int splicing_i_fraglen_check_cigar(const splicing_reads_t *reads, int idx) {
  const char *cigar=splicing_strvector_get(&reads->cigar, idx);
  char *ptr=(char *) cigar;
  long l=strtol(cigar, &ptr, 10L);
  return (ptr[0] == 'M' && ptr[1]=='\0') ? l : -1;
}

int splicing_estimate_fragment_length_file(const splicing_exonset_t *exons,
					   const char *readsfile,
					   splicing_vector_int_t *fraglen) {
  samfile_t *infile=0;
  bam1_t *read = bam_init1();
  char *mode="rb";
  char *myindexfile=0;
  bam_index_t *idx=0;
  FILE * ifp;
  int flen;
  int noexons=splicing_vector_int_size(&exons->seqid);
  char buf[500];
  const size_t buflen=sizeof(buf)/sizeof(char);
  int i, tid, beg, end, result;
  
  bam_verbose=0;
  
  SPLICING_FINALLY(splicing_i_bam_destroy1, read);
  
  splicing_vector_int_clear(fraglen);

  flen=strlen(readsfile);
  myindexfile = malloc(flen + 5);
  strcpy(myindexfile, readsfile);
  strcat(myindexfile, ".bai");
  ifp = fopen(myindexfile, "rb");
  free(myindexfile);
  if (ifp) { 
    idx = bam_index_load_core(ifp);
    SPLICING_FINALLY(bam_index_destroy, idx);
    fclose(ifp);
  }
  if (!idx) {
    SPLICING_ERROR("Cannot read BAM index file", SPLICING_EFILE);
  }

  infile=samopen(readsfile, mode, /*aux=*/ 0);
  if (!infile) { 
    SPLICING_ERROR("Cannot open SAM/BAM file", SPLICING_EFILE);
  }

  for (i=0; i<noexons; i++) {
    const char *seq=splicing_strvector_get(&exons->seqids, 
					   VECTOR(exons->seqid)[i]);
    snprintf(buf, buflen, "%s:%i-%i", seq, VECTOR(exons->start)[i],
	     VECTOR(exons->end)[i]);
    bam_parse_region(infile->header, buf, &tid, &beg, &end);
    if (tid < 0) {
      /* We silently skip this exon, it has no reads */
    } else {
      splicing_reads_t reads;
      int j, n, l, rl1, rl2;
      splicing_reads_init(&reads);
      result = bam_fetch(infile->x.bam, idx, tid, beg, end, &reads,
			 splicing_i_estfraglen_cb);
      if (reads.noPairs != 0) {
	SPLICING_CHECK(splicing_i_order_reads(&reads));	
      }
      n=reads.noPairs;
      for (j=0; j<n; ) {
	if (VECTOR(reads.pairpos)[j] <= 0) { 
	  j++; 
	} else {
	  int pos=VECTOR(reads.position)[j];
	  int ppos=VECTOR(reads.pairpos)[j];
	  /* Cigar strings must have a single matching part */
	  rl1=splicing_i_fraglen_check_cigar(&reads, j);
	  rl2=splicing_i_fraglen_check_cigar(&reads, j+1);
	  if (rl1 > 0 && rl2 > 0 && 
	      pos  >= beg && pos+rl1-1  <= end && 
	      ppos >= beg && ppos+rl2-1 <= end) {
	    l=ppos + rl2 - 1 - pos + 1;
	    SPLICING_CHECK(splicing_vector_int_push_back(fraglen, l));
	  }
	  j+=2;
	}
      }
      splicing_reads_destroy(&reads);
    }
  }
  
  return 0;
}

int splicing_estimate_fragment_length(const splicing_exonset_t *exons,
				      const char *readsfile,
				      const splicing_reads_t *reads,
				      splicing_vector_int_t *fraglen) {

  if ( (readsfile ? 1 : 0) + (reads ? 1 : 0) != 1) {
    SPLICING_ERROR("Give exactly one of `readsfile' and `reads'", 
		   SPLICING_EINVAL);
  }
  
  if (readsfile) { 
    return splicing_estimate_fragment_length_file(exons, readsfile, fraglen);
  } else {
    SPLICING_ERROR("`reads' is not yet implemented", 
		   SPLICING_UNIMPLEMENTED);
  }
}
