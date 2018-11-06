# Phylib: A library for phylogenetic tree inference

Author: David James Bryant david.bryant@otago.ac.nz

WE implement ant STL type container for iterator over and manipulating phylogenetic trees, together with a number of utitility and simulation algorithms. 

## Contents


*   `phylib.h`  Overall header file - include this file to use Phylib. 

*  global    Header files included in `phylib.h`
*  trees
   * `phylo.h`  Main STL-type class for phylogenies
   * `phylo_io.h` Input-output routines for the phylo class. 
   * `phylo_sim.h` Generation of random phylogenies from different distributions
   *  `phylo_util.h` General tools for manipulating and analysing phylogenies
   *  `menman.h` Memory management tools used by `phylo.h`
*  utilities   Generic useful utilities, not necessarily related to trees.



## Building

Include `phylo.h`and compile. **TODO**

## Overall approach


The phylogeny is seen as a container to hold information at the nodes, and the type of information stored is specified using the template. For example, a `phylo<int>` would store a single integer for every node. This makes it easy to write tree based algorithms: decide on what information is to be stored at a node and include this in the class. For example, an algorithm for Fitch parsimony might store:
```
	   class {
	   		BitSet states;
	   		int taxon id;
	   }	 
```
at the nodes. 
	   
 Unlike lists and vectors and so on, the structure of the container (i.e. the underlying tree) is important. We therefore add convertors to go from two different types of phylo: to copy from `phylo<type1>` to `phylo<type2>`. If you are implementing several algorithms one program you can  either define a different phylo type for each algorithm, or put the data from all the algorithms into one class. The advantage of the former is simplicty and memory reduction. The advantage of the latter is that you can reduce the number of timesphylogenies are copied.
	   
We note that the ordering of children is maintained, even though it (technically) doesn't affect the meaning of the phylogeny. Thus two
phylogenies are equal if and only if the trees are the same, their contents are the same, and the ordering of the children of 
each node is the same. See the tree utilities section for routines to compare trees that do not respect child ordering (a much more
computationally difficult issue).

The storage of the tree is based on a data structure designed (I think) by Dave Swofford, and used to great effect in PAUP. There is an
'immortal' header node that can never be removed. It usually has exactly one child, the root of the tree.[CHECK]
A node in the tree has three links:
 * one to its parent (or to the header if this is the root)
 * one to its leftmost child (or to null if this node has no children)
 * one to its next sibling (or to null if this is the rightmost child)
I have programmed using both this data structure and one where children are stored as an array, and find that Dave's structure
is easier to work with, after a bit of a steep learning curve. In future, I'll implement a binary tree only version that just
stores the left and right child.
	   
The basic phylo class provides a single iterator. Increments  `p++`are **not** implemented, as they are not uniquely defined. Instead we provide two ways to advance the iterator: `next_pre()`and `next_post`which implement a pre-order traversal (every node is visited before its children) and a post order traversal (each node is visited after its children). Both return a `null`iterator  (i.e, `p.null()==true`) if there are no further nodes in the traversal. There is alo an option to specify the node equal to the root of a subtree being traversed.



  
   I based the implementation of the phylo on the SGI list code, as adapted by D.Musser for a tutorial on lists. As such, I include the
   HP copyright header. Once I figure out how, I'll place this code under the GNU lesser license.
   
  
   Copyright (c) 1994
   Hewlett-Packard Company
  
   Permission to use, copy, modify, distribute and sell this software
   and its documentation for any purpose is hereby granted without fee,
   provided that the above copyright notice appear in all copies and
   that both that copyright notice and this permission notice appear
   in supporting documentation.  Hewlett-Packard Company makes no
   representations about the suitability of this software for any
   purpose.  It is provided "as is" without express or implied warranty.
  
   Adapted by D.Musser for a tutorial on lists.
  
   @version 0.2 7/11/18 Initial implementation by DJB
	  /



