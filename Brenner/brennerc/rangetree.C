/*
 This code was released into the public domain by Peter McCluskey on 8/24/2000.
*/

#include "brenner.h"
#include <stdlib.h>
#include <assert.h>

#ifndef __CCAPI_H__
#include "ccapi.h"
#endif

/* calc_dist needs State to be defined, so it has to appear after ccapi.h. */
#ifndef __CALC_DIST_H__
#include "calc_dist.h"
#endif

#ifndef __CCNEIGHBORSTATE_H__
#include "CCNeighborState.h"
#endif

/*
 * The following code is a way of producing O(N log^2 N) scaleup for large
 * molecules. Large in this context means it has atoms which are farther
 * apart than the rslj cutoff, which typically happens around 5000 to
 * 10000 atoms if the molecule is fairly compact. The execution time will
 * scale up at O(N^2) until that cutoff becomes important, and there isn't
 * much that can be done about that without sacrificing some quality in
 * the results.
 *
 * The approach used here is the bridged (layered) range-tree method
 * described on pages 83 through 88 of Preparata & Shamos' Computational
 * Geometry.
 */

typedef struct struct2dtree
{
  Double min_y;		/* if you're short on RAM, changing these to float */
  Double max_y;		/* may speed it up */
  struct struct2dtree *lower;
  struct struct2dtree *upper;
  const BrenAtom **atoms;
  int num_atoms;
  int *indices;
} Tree2d;

typedef struct struct3dtree
{
  IntPair index_ranges;
  Double min_x;
  Double max_x;
  struct struct3dtree *lower;
  struct struct3dtree *upper;
  Tree2d *tree2d;
  const BrenAtom **atoms;
} Tree3d;


static int cmp_number(const void *v1, const void *v2) /* probably temp hack */
{
  const BrenAtom *a1 = *(const BrenAtom**)v1;
  const BrenAtom *a2 = *(const BrenAtom**)v2;
  if(a1->number < a2->number) return -1;
  if(a1->number > a2->number) return 1;
  return 0;
}

static int cmp_x(const void *v1, const void *v2)
{
  const BrenAtom *a1 = *(const BrenAtom**)v1;
  const BrenAtom *a2 = *(const BrenAtom**)v2;
#if 1  /* uglier but faster */
  static const int c[3] = { -1, 1, 0 };
  return c[(a1->coord.x >= a2->coord.x) + (a1->coord.x == a2->coord.x)];
#else
  Double d = a1->coord.x - a2->coord.x;
  if(d < 0) return -1;
  if(d > 0) return 1;
  return 0;
#endif
}

static int cmp_y(const void *v1, const void *v2)
{
  const BrenAtom *a1 = *(const BrenAtom**)v1;
  const BrenAtom *a2 = *(const BrenAtom**)v2;
#if 1  /* uglier but faster */
  static const int c[3] = { -1, 1, 0 };
  return c[(a1->coord.y >= a2->coord.y) + (a1->coord.y == a2->coord.y)];
#else
  Double d = a1->coord.y - a2->coord.y;
  if(d < 0) return -1;
  if(d > 0) return 1;
  return 0;
#endif
}

static int cmp_z(const void *v1, const void *v2)
{
  const BrenAtom *a1 = *(const BrenAtom**)v1;
  const BrenAtom *a2 = *(const BrenAtom**)v2;
#if 1  /* uglier but faster */
  static const int c[3] = { -1, 1, 0 };
  return c[(a1->coord.z >= a2->coord.z) + (a1->coord.z == a2->coord.z)];
#else
  Double d = a1->coord.z - a2->coord.z;
  if(d < 0) return -1;
  if(d > 0) return 1;
  return 0;
#endif
}

static int
squeeze_atoms_x(const BrenAtom **atoms, int num_atms, Double max_cutoff,
                const BrenAtom **result)
{
  int i, j;
  int k = 0;
  int last_movable = -1;
  for(i = 0; i < num_atms; ++i)
    {
      if(atoms[i]->movable || i == num_atms-1)
        {
          for(j = last_movable + 1; j < i; ++j)
            {
              if(fabs(atoms[j]->coord.x - atoms[i]->coord.x) <= max_cutoff
                 || (last_movable != -1
                     && fabs(atoms[j]->coord.x - atoms[last_movable]->coord.x) <= max_cutoff))
                {
                  result[k++] = atoms[j];
                }
            }
          if(atoms[i]->movable)
            {
              last_movable = i;
              result[k++] = atoms[i];
            }
        }
    }
  return k;
}

static int
squeeze_atoms_y(const BrenAtom **atoms, int num_atms, Double max_cutoff,
                const BrenAtom **result)
{
  int i, j;
  int k = 0;
  int last_movable = -1;
  for(i = 0; i < num_atms; ++i)
    {
      if(atoms[i]->movable || i == num_atms-1)
        {
          for(j = last_movable + 1; j < i; ++j)
            {
              if(fabs(atoms[j]->coord.y - atoms[i]->coord.y) <= max_cutoff
                 || (last_movable != -1
                     && fabs(atoms[j]->coord.y - atoms[last_movable]->coord.y) <= max_cutoff))
                {
                  result[k++] = atoms[j];
                }
            }
          if(atoms[i]->movable)
            {
              last_movable = i;
              result[k++] = atoms[i];
            }
        }
    }
  return k;
}

static int
squeeze_atoms_z(const BrenAtom **atoms, int num_atms, Double max_cutoff,
                const BrenAtom **result)
{
  int i, j;
  int k = 0;
  int last_movable = -1;
  for(i = 0; i < num_atms; ++i)
    {
      if(atoms[i]->movable || i == num_atms-1)
        {
          for(j = last_movable + 1; j < i; ++j)
            {
              if(fabs(atoms[j]->coord.z - atoms[i]->coord.z) <= max_cutoff
                 || (last_movable != -1
                     && fabs(atoms[j]->coord.z - atoms[last_movable]->coord.z) <= max_cutoff))
                {
                  result[k++] = atoms[j];
                }
            }
          if(atoms[i]->movable)
            {
              last_movable = i;
              result[k++] = atoms[i];
            }
        }
    }
  return k;
}

static int log2int(int n)
{
  int lg2 = 0;
  assert(n >= 2);
  --n;
  do ++lg2; while(n >>= 1);
  return lg2;
}

int cnt_alloc[3];

static void *
build2dtree(Tree2d *tree2, const BrenAtom **atom_list,
            const int *atom_number, int depth, IntPair index_ranges,
            const BrenAtom **zrange, char *mem2d_block)
{
  int i;
  int diff = index_ranges.v2 - index_ranges.v1;
  tree2->num_atoms = diff;
  tree2->min_y = atom_list[index_ranges.v1]->coord.y;
  tree2->max_y = atom_list[index_ranges.v1 + diff-1]->coord.y;
  if(diff > 1)
    {
      const BrenAtom **child_atoms[2];
      IntPair lower_range, upper_range;
      int median = index_ranges.v1 + diff/2;
      int n1 = diff/2;
      int n2 = index_ranges.v2 - median;
      int *lower_indices, *upper_indices;
      int child_index[2];
      child_index[0] = child_index[1] = 0;
      tree2->atoms = zrange;
      tree2->lower = (Tree2d *)mem2d_block;
      tree2->upper = tree2->lower + 1;
      mem2d_block += 2*sizeof(Tree2d);
      cnt_alloc[0] += 2*sizeof(Tree2d);
      
      child_atoms[0] = (const BrenAtom **)mem2d_block;
      mem2d_block += n1*sizeof(const BrenAtom*);
      child_atoms[1] = (const BrenAtom **)mem2d_block;
      mem2d_block += n2*sizeof(const BrenAtom*);
      cnt_alloc[1] += n2*sizeof(const BrenAtom*) + n1*sizeof(const BrenAtom*);
      tree2->lower->indices = (int *)mem2d_block;
      mem2d_block += (diff+1)*sizeof(tree2->indices[0]);
      tree2->upper->indices = (int *)mem2d_block;
      mem2d_block += (diff+1)*sizeof(tree2->indices[0]);
      cnt_alloc[2] += 2*(diff+1)*sizeof(tree2->indices[0]);
#if 0
      if(0)
        printf("%*salloc%*s %6d %8d %5d\n", 
               depth, "", 9-depth, "", mem2d_block - save_mem2d_block,
               mem2d_block - mem_orig, diff);
      if(mem2d_block >= mem_orig + block2d_size)
        {
          int i;
          static int sum[3];
          for(i = 1; i < 20; ++i)
            {
              printf("%2d %8d %8d %8d\n", i, cnt_alloc[0][i], cnt_alloc[1][i], cnt_alloc[2][i]);
              sum[0] += cnt_alloc[0][i];
              sum[1] += cnt_alloc[1][i];	
              sum[2] += cnt_alloc[2][i];
            }
          printf("\t%d %d %d  %d\n", sum[0], sum[1], sum[2], sum[0] + sum[1] + sum[2]);
          abort();
        }
      fflush(stdout);
#endif
      
      lower_indices = tree2->lower->indices;
      upper_indices = tree2->upper->indices;
      for(i = 0; i < diff; ++i)
        {
          const BrenAtom *p = tree2->atoms[i];
          int childno = (atom_number[p->number] >= median);
          lower_indices[i] = child_index[0];
          upper_indices[i] = child_index[1];
          child_atoms[childno][child_index[childno]++] = p;
        }
      lower_indices[i] = child_index[0];
      upper_indices[i] = child_index[1];
      assert(n1 == child_index[0] && n2 == child_index[1]);
      lower_range.v1 = index_ranges.v1;
      lower_range.v2 = median;
      mem2d_block = build2dtree(tree2->lower, atom_list, atom_number, depth + 1,
                                lower_range, child_atoms[0], mem2d_block);
      upper_range.v1 = median;
      upper_range.v2 = index_ranges.v2;
      mem2d_block = build2dtree(tree2->upper, atom_list, atom_number, depth + 1,
                                upper_range, child_atoms[1], mem2d_block);
    }
  else			/* at leaf */
    {
      tree2->atoms = zrange;
      tree2->lower = NULL;
      tree2->upper = NULL;
    }
  return mem2d_block;
}

#if 0
static void verify2(const Tree2d *tree2)
{
  if(tree2->lower)
    verify2(tree2->lower);
  if(tree2->upper)
    verify2(tree2->upper);
}

static void verify3(const Tree3d *tree3)
{
  if(tree3->tree2d)
    verify2(tree3->tree2d);
  if(tree3->lower)
    verify3(tree3->lower);
  if(tree3->upper)
    verify3(tree3->upper);
}
#endif

static void
build3dtree(Tree3d *tree3, const BrenAtom **atom_list,
	    const BrenAtom **yrange, const BrenAtom **zrange,
	    const int *atom_xnum, int depth, int max_atoms)
{
  int i;
  int diff = tree3->index_ranges.v2 - tree3->index_ranges.v1;
  tree3->min_x = atom_list[tree3->index_ranges.v1]->coord.x;
  assert(tree3->index_ranges.v2 > 0);
  tree3->max_x = atom_list[tree3->index_ranges.v2-1]->coord.x;
  if(diff > 1)
  {
    int median = tree3->index_ranges.v1 + diff/2;
    IntPair index_range2;
    static int *atom_number;
    static int *atom_znumber;

    static char *mem2d_block_static;
    static int mem2d_block_index;
    static int mem2d_block_allocated;

    int log2n = log2int(diff);
    int block2d_size = 2*diff*sizeof(Tree2d)
      + (log2n + 1)*diff*sizeof(const BrenAtom*)
      + (2*log2n + 2)*diff*sizeof(tree3->tree2d->indices[0]);
    char *mem2d_block;

    const BrenAtom **lower_atoms = NULL;
    const BrenAtom **upper_atoms = NULL;
    const BrenAtom **lower_zrange = NULL;
    const BrenAtom **upper_zrange = NULL;
    static Tree3d *tree_array;
    static int tree_array_index;
    static char *mem3d_block;
    static int mem3d_block_index;
    static int mem3d_block_allocated;
    char *mem2d_block_limit;
    if(!depth)
    {
      static int last_num_atoms;
      static int last_max_atoms;
      if(max_atoms != last_max_atoms)
      {
	last_max_atoms = max_atoms;
	if(atom_number) free(atom_number);
	atom_number = (int *)malloc(max_atoms * sizeof(atom_number[0]));
	if(atom_znumber) free(atom_znumber);
	atom_znumber = (int *)malloc(max_atoms * sizeof(atom_znumber[0]));
      }
      tree_array_index = 0;
      mem3d_block_index = 0;
      mem2d_block_index = 0;
      if(diff > last_num_atoms)
      {
	if(mem3d_block) free(mem3d_block);
	if(tree_array) free(tree_array);
	tree_array = (Tree3d *)malloc(2*diff*sizeof(Tree3d));
	mem3d_block_allocated = 2*log2n*diff*sizeof(const BrenAtom*);
	mem3d_block = malloc(mem3d_block_allocated);
	if(mem2d_block_static) free(mem2d_block_static);
	mem2d_block_allocated = log2n * block2d_size;
	mem2d_block_static = malloc(mem2d_block_allocated);
#ifndef NDEBUG
	printf("%d alloced for build3dtree, %d for 2d (%d*%d)\n",
	       mem3d_block_allocated, mem2d_block_allocated, log2n, block2d_size);
#endif
	last_num_atoms = diff;
      }
    }
    mem2d_block = mem2d_block_static + mem2d_block_index;
    mem2d_block_index += block2d_size;
    mem2d_block_limit = mem2d_block + block2d_size;
    assert(mem2d_block_index <= mem2d_block_allocated);
    if(0)
      printf("log2n %3d block2d_size %5d %5d %5d %5d\n",
	     log2n, block2d_size, 2*diff*sizeof(Tree2d),
	     (log2n + 1)*diff*sizeof(const BrenAtom*),
	     (2*log2n + 2)*diff*sizeof(tree3->tree2d->indices[0]));

    assert(yrange);
    tree3->atoms = yrange;
    {
      int n1 = median - tree3->index_ranges.v1;
      int n2 = tree3->index_ranges.v2 - median;
      int i1 = 0, i2 = 0, z1 = 0, z2 = 0;
      if(n1 > 1)
      {
	lower_atoms = (const BrenAtom **)(mem3d_block + mem3d_block_index);
	mem3d_block_index += n1*sizeof(const BrenAtom*);
	lower_zrange = (const BrenAtom **)(mem3d_block + mem3d_block_index);
	mem3d_block_index += n1*sizeof(const BrenAtom*);
      }
      else lower_atoms = lower_zrange = &atom_list[tree3->index_ranges.v1];
      if(n2 > 1)
      {
	upper_atoms = (const BrenAtom **)(mem3d_block + mem3d_block_index);
	mem3d_block_index += n2*sizeof(const BrenAtom*);
	upper_zrange = (const BrenAtom **)(mem3d_block + mem3d_block_index);
	mem3d_block_index += n2*sizeof(const BrenAtom*);
      }
      else upper_atoms = upper_zrange = &atom_list[median];
      assert(mem3d_block_index <= mem3d_block_allocated);
      for(i = 0; i < diff; ++i)
      {
	const BrenAtom *p = tree3->atoms[i];
	if(atom_xnum[p->number] < median)
	  lower_atoms[i1++] = p;
	else
	  upper_atoms[i2++] = p;
	p = zrange[i];
	if(atom_xnum[p->number] < median)
	  lower_zrange[z1++] = p;
	else
	  upper_zrange[z2++] = p;
      }
      assert(n1 == i1 && n2 == i2);
    }

    tree3->lower = &tree_array[tree_array_index++];
    tree3->lower->index_ranges.v1 = tree3->index_ranges.v1;
    tree3->lower->index_ranges.v2 = median;
    build3dtree(tree3->lower, atom_list, lower_atoms, lower_zrange,
		atom_xnum, depth + 1, max_atoms);
    tree3->upper = &tree_array[tree_array_index++];
    tree3->upper->index_ranges.v1 = median;
    tree3->upper->index_ranges.v2 = tree3->index_ranges.v2;
    build3dtree(tree3->upper, atom_list, upper_atoms, upper_zrange,
		atom_xnum, depth + 1, max_atoms);

    for(i = 0; i < diff; ++i)
    {
      assert(yrange[i]->number < max_atoms);
      assert(zrange[i]->number < max_atoms);
      atom_number[yrange[i]->number] = i;
      atom_znumber[zrange[i]->number] = i;
    }
    index_range2.v1 = 0;
    index_range2.v2 = diff;
    tree3->tree2d = (Tree2d *)mem2d_block;
    mem2d_block += sizeof(Tree2d);
    tree3->tree2d->indices = NULL;
    mem2d_block = build2dtree(tree3->tree2d, tree3->atoms, atom_number,
			      0, index_range2, zrange, mem2d_block);
#if 0
      printf("%6d of %7d unused %5d atoms %5d %5d %5d %p %p %6d\n",
	     mem2d_block_limit - mem2d_block,
	     block2d_size, diff, cnt_alloc[0], cnt_alloc[1], cnt_alloc[2],
	     mem2d_block, mem2d_block_limit, mem2d_block_index);
#if 0
    if(mem2d_block > mem2d_block_limit)
      for(i = 0; i < 11; ++i) printf("%3d %5d\n", i, cnt22[i]);
#endif
#endif
    assert(mem2d_block <= mem2d_block_limit);
  }
  else			/* at leaf */
  {
    tree3->atoms = &atom_list[tree3->index_ranges.v1];
    if(0)
    printf("%*s (%2d-%2d) (%4.1f,%4.1f,%4.1f)\n",
	   depth,"", tree3->index_ranges.v1, tree3->index_ranges.v2,
	   atom_list[0]->coord.x,
	   atom_list[0]->coord.y,
	   atom_list[0]->coord.z);
    tree3->lower = NULL;
    tree3->upper = NULL;
    tree3->tree2d = NULL;
  }
}

static int
findz(const BrenAtom **atom_range, int num_atms, Double val, int dir)
{
  int i = num_atms/2;
  int delta = (num_atms+1)/2;
  while(delta > 1)
  {
    Double diff = atom_range[i]->coord.z - val;
    if(diff == 0) return i;
    delta = (delta + 1)/2;
    i += (diff < 0 ? 1 : -1) * delta;
    if(i < 0) i = 0;
    else if(i >= num_atms) i = num_atms - 1;
  }
  return i;
}

static int
search2d(const Tree2d *tree2, const BrenAtom **candidates,
	 int start_index, int stop_index,
	 const double min_coord[3], const double max_coord[3],
	 int min_atom_index, int depth)
{
  int sum = 0;
  int i;
  const Tree2d *lower = tree2->lower;
  if(lower)
  {
    int new_start, new_stop;
    const Tree2d *upper;
    if(depth == 0) /* called from search2droot */
    {
      new_start = start_index;
      new_stop  = stop_index;
    }
    else if(start_index == 0)
    {
      new_start = 0;
      new_stop = tree2->indices[stop_index];
    }
    else
    {
      new_start = tree2->indices[start_index];
      new_stop = tree2->num_atoms;
    }
    /* assert(new_stop <= tree2->upper->num_atoms); */
    if(new_start >= new_stop)
      return 0;

    upper = tree2->upper;
    if(min_coord[1] <= lower->min_y && max_coord[1] >= upper->max_y)
    {
      int n = tree2->num_atoms;
      /* assert(stop_index <= n); */
      for(i = new_start; i < new_stop; ++i)
	if(tree2->atoms[i]->number > min_atom_index)
	  candidates[sum++] = tree2->atoms[i];
      if(0 && new_stop > 0)
      printf("%*s found range (%2d[%4.1f] to %2d[%4.1f])  %2d  %2d\n",
	     depth, "", tree2->atoms[new_start]->number, tree2->atoms[new_start]->coord.z,
	     tree2->atoms[new_stop-1]->number, tree2->atoms[new_stop-1]->coord.z, sum, n);
      return sum;
    }
    if(min_coord[1] < upper->min_y && lower->min_y <= max_coord[1])
    {
      sum = search2d(lower, candidates, new_start, new_stop,
		     min_coord, max_coord, min_atom_index, depth + 1);
    }
    if(max_coord[1] > lower->max_y && upper->max_y >= min_coord[1])
    {
      sum += search2d(upper, &candidates[sum], new_start, new_stop,
		      min_coord, max_coord, min_atom_index, depth + 1);
    }
    if(0)
      printf("%*ssearch2d%*s%3d %3d %7.3f %7.3f  %7.3f %7.3f tot %2d %p\n",
	     depth, "", 9-depth, "", start_index, stop_index,
	     min_coord[1], tree2->lower->min_y,
	     max_coord[1], tree2->upper->max_y, sum, tree2);
    return sum;
  }
  else
  {
    if(tree2->num_atoms > 0)
    {
      const BrenAtom *p = tree2->atoms[0];
      if(p->number > min_atom_index
	 && p->coord.z >= min_coord[2] && p->coord.z <= max_coord[2])
      {
	candidates[0] = p;
	if(0) printf("candidate0 %3d\n", p->number);
	return 1;
      }
    }
  }
  return 0;
}

static int
search2droot(const Tree2d *tree2, const BrenAtom **candidates,
	     const double min_coord[3], const double max_coord[3],
	     int min_atom_index)
{
  int sum = 0;
  int i;
  const Tree2d *lower = tree2->lower;
  if(lower)
  {
    int n = tree2->num_atoms;
    int med_index;
    int is_low_end;
    if(min_coord[1] <= lower->min_y && max_coord[1] >= tree2->upper->max_y)
    {
      int start_index = findz(tree2->atoms, n, min_coord[2], -1);
      int stop_index  = findz(tree2->atoms, n, max_coord[2], 1);
      if(stop_index < n && tree2->atoms[stop_index]->coord.z < max_coord[2])
	++stop_index;
      assert(stop_index <= n);
      if(start_index && tree2->atoms[start_index]->coord.z > min_coord[2])
	--start_index;
      for(i = start_index; i < stop_index; ++i)
	if(tree2->atoms[i]->number > min_atom_index)
	  candidates[sum++] = tree2->atoms[i];
      if(0)
	printf(" root found range (%2d[%4.1f] to %2d[%4.1f])  %2d  %2d\n",
	     tree2->atoms[start_index]->number, tree2->atoms[start_index]->coord.z,
	     tree2->atoms[stop_index-1]->number, tree2->atoms[stop_index-1]->coord.z, sum, n);
      return sum;
    }
    med_index = -1;
    is_low_end = (min_coord[2] <= tree2->atoms[0]->coord.z);
    if(min_coord[1] < tree2->upper->min_y && lower->min_y <= max_coord[1])
    {
      int new_start, new_stop;
      if(is_low_end)
      {
	new_start = 0;
	med_index = findz(tree2->atoms, n, max_coord[2], 1);
	if(med_index < n && tree2->atoms[med_index]->coord.z < max_coord[2])
	  ++med_index;
	new_stop = lower->indices[med_index];
      }
      else
      {
	med_index = findz(tree2->atoms, n, min_coord[2], -1);
	if(med_index && tree2->atoms[med_index]->coord.z > min_coord[2])
	  --med_index;
	new_start = lower->indices[med_index];
	new_stop = lower->num_atoms;
      }
      if(new_start < new_stop)
	sum = search2d(lower, candidates, new_start, new_stop,
		       min_coord, max_coord, min_atom_index, 0);
    }
    if(max_coord[1] > lower->max_y && tree2->upper->max_y >= min_coord[1])
    {
      int new_start, new_stop;
      if(is_low_end)
      {
	new_start = 0;
	if(med_index == -1)
	{
	  med_index = findz(tree2->atoms, n, max_coord[2], 1) + 1;
	  if(med_index < n && tree2->atoms[med_index]->coord.z < max_coord[2])
	    ++med_index;
	}
	new_stop = tree2->upper->indices[med_index];
      }
      else
      {
	if(med_index == -1)
	{
	  med_index = findz(tree2->atoms, n, min_coord[2], -1);
	  if(med_index && tree2->atoms[med_index]->coord.z > min_coord[2])
	    --med_index;
	}
	new_start = tree2->upper->indices[med_index];
	new_stop = tree2->upper->num_atoms;
      }
      if(new_start < new_stop)
	sum += search2d(tree2->upper, &candidates[sum], new_start, new_stop,
			min_coord, max_coord, min_atom_index, 0);
    }
    if(0)
      printf("search2d %3d %3d %7.3f %7.3f   %7.3f %7.3f tot %2d Root %p\n",
	     med_index, n,
	     min_coord[1], lower->min_y,
	     max_coord[1], tree2->upper->max_y, sum, tree2);
    return sum;
  }
  else
  {
    if(tree2->num_atoms > 0)
    {
      const BrenAtom *p = tree2->atoms[0];
      if(p->number > min_atom_index
	 && p->coord.z >= min_coord[2] && p->coord.z <= max_coord[2])
      {
	candidates[0] = p;
	if(0) printf("candidate0 %3d root\n", p->number);
	return 1;
      }
    }
  }
  return 0;
}

static int
search3d(const Tree3d *tree3, const BrenAtom **candidates,
	 const double min_coord[3], const double max_coord[3],
	 int min_atom_index, int depth)
{
  int sum = 0;
  int i;
  if(tree3->lower)
  {
    if(min_coord[0] <= tree3->lower->min_x && max_coord[0] >= tree3->upper->max_x)
    {
      if(0)
	printf("%*sdescend%*sto 2d %2d %2d (%4.1f %4.1f) %4.1f %4.1f\n",
	     depth, "", 12-depth, "",
	     tree3->index_ranges.v1, tree3->index_ranges.v2,
	     tree3->lower->min_x, tree3->upper->max_x, min_coord[0], max_coord[0]);
      if(tree3->tree2d)
	return search2droot(tree3->tree2d, candidates,
			    min_coord, max_coord, min_atom_index);
    }
    if(min_coord[0] < tree3->upper->min_x)
      sum = search3d(tree3->lower, candidates, min_coord, max_coord,
		      min_atom_index, depth + 1);
    if(max_coord[0] > tree3->lower->max_x)
      sum += search3d(tree3->upper, &candidates[sum], min_coord, max_coord,
		      min_atom_index, depth + 1);
    if(0)
      printf("%*sgot%*s%2d %7.3f %7.3f   %7.3f %7.3f (range %6.3f %6.3f) %p\n",
	     depth, "", 12-depth, "", sum, min_coord[0], tree3->upper->min_x,
	     max_coord[0], tree3->lower->max_x, tree3->min_x, tree3->max_x, tree3);
  }
  else
  {
    const BrenAtom *a = tree3->atoms[0];
    if(0)
      printf("%*s3leaf %d %4.1f (%4.1f - %4.1f) %4.1f (%4.1f - %4.1f)\n",
	   depth, "", a->number, a->coord.y, min_coord[1], max_coord[1],
	   a->coord.z, min_coord[2], max_coord[2]);
    if(a->number > min_atom_index
       && a->coord.y >= min_coord[1] && a->coord.y <= max_coord[1]
       && a->coord.z >= min_coord[2] && a->coord.z <= max_coord[2])
    {
      candidates[0] = a;
      return 1;
    }
  }
  if(0)
  printf("%*ssearch3d%*s(%4.1f - %4.1f) got %d (%2d,%2d)\n", depth,"",
	 12-depth, "",tree3->min_x,tree3->max_x, sum,
	 tree3->index_ranges.v1,tree3->index_ranges.v2);
  return sum;
}

/* Compute neighbors and put the result into ljns.  ljns may be the
   lennard-jones neighbor structure or the bond-order neighbor structure.  rslj
   says how far apart atoms can be and still be neighbors, as a function of the
   ktype's of the two atoms.  xmms is the minimum cutoff; if two atoms are
   closer than xmms then we don't want to consider them to be neighbors
   either. */
void
make_neighbor_lists_tree(BrennerMainInfo *s, struct NeighborState *ljns,
			 Float rslj[NTYPES+1][NTYPES+1],
			 Float xmms[NTYPES+1][NTYPES+1])
{
  int i,j,k,n,x,y,z;
  int index1;
  int kend = 0;
	/* num_atoms_used is how many atoms are within a distance of
	   max_cutoff from a movable atom along x, y, or z axis */
  int num_atoms_used;
  const int num_atms = s->num_atms;
  const BrenAtom *const atm_num = s->atm_num;
  Tree3d tree3;
  static const BrenAtom **xrange;
  static const BrenAtom **yrange;
  static const BrenAtom **zrange;
  static int *atom_xnum;
    /* candidates will hold pointers to the neighbors themselves. */
  static const BrenAtom **candidates;
  static int *last_neighbor_found;
  static int last_num_atoms;
  Double max_cutoff = 0;
  allocate_atom (ljns, s->num_atms);
  if(num_atms > last_num_atoms)
  {
    last_num_atoms = num_atms;
    if(xrange) free(xrange);
    if(yrange) free(yrange);
    if(zrange) free(zrange);
    if(atom_xnum) free(atom_xnum);
    if(candidates) free(candidates);
    if(last_neighbor_found) free(last_neighbor_found);
    xrange = (const BrenAtom **)malloc(num_atms * sizeof(*xrange));
    yrange = (const BrenAtom **)malloc(num_atms * sizeof(*yrange));
    zrange = (const BrenAtom **)malloc(num_atms * sizeof(*zrange));
    atom_xnum = (int *)malloc(num_atms * sizeof(*atom_xnum));
    candidates = (const BrenAtom **)malloc(num_atms * sizeof(*candidates));
    last_neighbor_found = (int *)malloc(num_atms * sizeof(*last_neighbor_found));
  }
  for(i = 1; i <= NTYPES; ++i)
  {
    for(j = 1; j <= NTYPES; ++j)
    {
        if(rslj[i][j] > max_cutoff) max_cutoff = rslj[i][j];
    }
  }
  for(i = 0; i < num_atms; ++i)
  {
    xrange[i] = &atm_num[i];
    last_neighbor_found[i] = -1;
  }

  if(s->squeeze)
  {
    qsort(xrange, num_atms, sizeof(xrange[0]), cmp_z);
    num_atoms_used = squeeze_atoms_z(xrange, num_atms, max_cutoff, zrange);
    qsort(zrange, num_atoms_used, sizeof(zrange[0]), cmp_y);
    num_atoms_used = squeeze_atoms_y(zrange, num_atoms_used, max_cutoff, yrange);
    qsort(yrange, num_atoms_used, sizeof(yrange[0]), cmp_x);
    num_atoms_used = squeeze_atoms_x(yrange, num_atoms_used, max_cutoff, xrange);
  }
  else
  {
    num_atoms_used = num_atms; 
    qsort(xrange, num_atoms_used, sizeof(xrange[0]), cmp_x);
  }

  for(i = 0; i < num_atms; ++i)
    atom_xnum[i] = -1;
  for(i = 0; i < num_atoms_used; ++i)
  {
    atom_xnum[xrange[i]->number] = i;
    yrange[i] = xrange[i];
    zrange[i] = xrange[i];
  }
  qsort(zrange, num_atoms_used, sizeof(zrange[0]), cmp_z);
  qsort(yrange, num_atoms_used, sizeof(yrange[0]), cmp_y);

  if(0) printf("num_atms %d num_atoms_used %d\n", num_atms, num_atoms_used);
  tree3.index_ranges.v1 = 0;
  tree3.index_ranges.v2 = num_atoms_used;
  build3dtree(&tree3, xrange, yrange, zrange, atom_xnum, 0, s->num_atms);

  qsort(xrange, num_atoms_used, sizeof(xrange[0]), cmp_number);/* FIXME must be better way */
  for(i = 0; i < num_atms; ++i)	{
    const BrenAtom *pAtm = &atm_num[i];
    const int ki = pAtm->ktype;
    const Float *rslj_ki = rslj[ki];
    const Float *xmms_ki = xmms[ki];
    /* n will be the number of potential neighbors found for atom i by
       searching the range tree. */
    int n;
    double min_coord[3];
    double max_coord[3];
    double max_r = 0;
    int j;
    ljns->neighbor_start [i] = kend;
    if(atom_xnum[i] == -1)
      continue;
    /* This next loop sets max_r to the maximum interaction distance between
       atom i and other atoms.  NTYPES is quite small, so we could have
       precomputed this. */
    for(j = 1; j <= NTYPES; ++j) {
      if(rslj[ki][j] > max_r) max_r = rslj[ki][j];
      if(!j) /* don't remove! without it, gcc version 2.95.2 with -O3 produces
                NaN for max_r under conditions I haven't fully tracked down;
                probably compiler bug */ 
        printf("rslj[%d][%d] %f\n",ki, j,  max_r);
    }
#if 0
    max_r = 0.0001; printf("hack %f\n", max_r);
#endif
    min_coord[0] = pAtm->coord.x - max_r;
    max_coord[0] = pAtm->coord.x + max_r;
    min_coord[1] = pAtm->coord.y - max_r;
    max_coord[1] = pAtm->coord.y + max_r;
    min_coord[2] = pAtm->coord.z - max_r;
    max_coord[2] = pAtm->coord.z + max_r;
    n = search3d(&tree3, candidates, min_coord, max_coord,
		 ljns == s->ljNeighborState ? pAtm->number : -1, 0);
    if(0)
      printf("search3d %d %d (%4.1f,%4.1f,%4.1f) max_r %6.2f\n",pAtm->number, n,
             pAtm->coord.x, pAtm->coord.y, pAtm->coord.z, max_r);
#if 0
    if(n != 1 || candidates[0]->number != i)
      exit(0);
#endif
    allocate_neighbor (ljns, kend + n);
    for(j = 0; j < n; ++j) {
      const BrenAtom *pAtm2 = candidates[j];
      int index2;
      int kj = pAtm2->ktype;
      double rsqs;
      if(rslj_ki[kj] == 0.0) continue;
      rsqs = CALC_DIST_NO_RR(toState (s), i, pAtm2->number);
      if(-1 == pAtm->number)
	printf("cmp %3d %3d %3d %d %d %7.2f %7.2f %5.1f %5.1f\n",i,
	       pAtm->number, pAtm2->number, ki, kj, rsqs, xmms_ki[kj],
	       posy (i, toState(s)), posy (pAtm2->number, toState(s)));
      if(rsqs > rslj_ki[kj]) continue;
      if(rsqs < xmms_ki[kj]) continue;
      index2 = pAtm2->number;
      if (index2 != pAtm->number && last_neighbor_found[index2] != pAtm->number) { 
	assert (kend < ljns->neighbors_allocated);
        ljns->neighbor_list [kend] = index2;
        ++kend;
	last_neighbor_found[index2] = pAtm->number;
      }
    }
  }
  assert (i == num_atms);
  ljns->neighbor_start [i] = kend;
  // There used to be code right here to sort the neighbor list.  But now the
  // neighbor list is organized differently, so that code doesn't work any
  // more.  We could sort the neighbors of each atom.  The only use for sorting
  // that I can see is debugging, so let's recreate the sorting later if and
  // when it's needed.  Tim Freeman 23 Aug 2000
}
