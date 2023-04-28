/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  tMeshList.h: Header file for derived classes tMeshList and tMeshListIter
**
**  A tMeshList is derived from the generic linked list class tList.
**  It is used in CHILD to store lists of grid elements (nodes and edges),
**  and differs from a generic list in being divided into two parts:
**  (1) an "active" part, representing elements that are not part of the
**  mesh boundary and are therefore subject to active processes (whatever
**  those may be; in CHILD the processes are runoff, erosion, and
**  sedimentation); and (2) a "boundary" part, containing elements along
**  the mesh boundary.
**
**  A tMeshListIter is an iterator for a tMeshList. It has the same services
**  as a tListIter. It also will move to the last "active" (non-boundary) 
**  node on a grid list, or to the first boundary node on the list. It adds 
**  special functions FirstP, NextP that are identical to the tListIter 
**  functions First and Next except that they return a pointer to the data
**  portion of the node (or zero if the end of the list is reached, or the 
**  current node is null for some other reason).
**
**************************************************************************/

#ifndef TMESHLIST_H
#define TMESHLIST_H

#include "src/Headers/Classes.h"
#include "src/Headers/Definitions.h"
#include "src/tList/tList.h"

//=========================================================================
//
//
//                  Section 1: tMeshList Class Declarations
//
//
//=========================================================================

/**************************************************************************
**
** tMeshList()
**
** Class tMeshList implements a linked list that is divided into two
** parts, an "active" (front) and "inactive" (back) part. Derived from tList.
**
**************************************************************************/

template< class NodeType >
class tMeshList 
	: public tList< NodeType >
{
  friend class tListIter< NodeType  >;
  friend class tMeshListIter< NodeType  >;

public:
  tMeshList();
  tMeshList( const tMeshList< NodeType > * );
  ~tMeshList();
  const tMeshList< NodeType >
       &operator=( const tMeshList< NodeType > & );
  int operator==( const tMeshList< NodeType > & ) const;
  int operator!=( const tMeshList< NodeType > & ) const;
  int getActiveSize() const
  { return nActiveNodes; }

#ifdef PARALLEL_TRIBS
  int getGlobalActiveSize();
  int* collectActiveSize();
#endif


  tListNode< NodeType  > * getLastActive() const;
  int isActiveEmpty() const;
  int isBoundEmpty() const;
  void insertAtBoundFront( const NodeType & );
  int removeFromBoundFront( NodeType & );
  void insertAtActiveBack( const NodeType & );
  int removeFromActiveBack( NodeType & );
  void setNActiveNodes( int );
  int removeNext( NodeType &value, tListNode< NodeType > * );
  int removePrev( NodeType &value, tListNode< NodeType > * );
  void moveToBack( tListNode< NodeType > * );
  void moveToFront( tListNode< NodeType > * );
  void moveToActiveBack( tListNode< NodeType > * );
  void moveToBoundFront( tListNode< NodeType > * );
  void moveToBack( NodeType * );
  void insertAtFront( const NodeType & );
  int removeFromFront( NodeType & );
  int InActiveList( tListNode< NodeType > * );
  void Flush();

  //SMM - added 08132008
  // Move next tListNode to back of list
  int nextToBack(tListNode<NodeType>*);
  // Move first tListNode to back of list
  int frontToBack();
   
protected:
  int nActiveNodes;                    // # of active nodes on list
  tListNode< NodeType > * lastactive;  // ptr to last active node
};

/**************************************************************************
** 
** tMeshListIter()
**
** Helper class for tMeshList, derived from tListIter ("iterators" that
** walk up and down a tList, fetching items). In addition to tListIter 
** capabilities, tMeshListIter adds methods to move to and/or fetch the 
** last active or first boundary (inactive) items, and to indicate whether 
** it is on currently on the active portion of the list.
**
**************************************************************************/

template< class NodeType >
class tMeshListIter : public tListIter< NodeType >
{
public: 
  tMeshListIter();
  tMeshListIter( tMeshList< NodeType > & );
  tMeshListIter( tMeshList< NodeType > *ptr )
        : tListIter< NodeType >( ptr )
  {
    assert( ptr != 0 );
    this->curnode = this->listPtr->first;
    assert( this->curnode != 0 );
  }
  ~tMeshListIter() {}

  int LastActive();
  int FirstBoundary();
  int IsActive()
  {
    int act;
    if( this->curnode!=0 ){
      assert( this->curnode->getDataPtr()!=0 );
      act = this->curnode->getDataRef().getBoundaryFlag();			
      if((act == kNonBoundary) || (act == kStream )) return 1;
    }
    return 0;
  }


  NodeType * LastActiveP();
  NodeType * FirstBoundaryP();
};

template< class NodeType >
int tMeshList< NodeType >::
removeNext( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr->next == 0 ) return 0;
   if( ptr == 0 ) return 0;
   if( ptr->next == lastactive ) return removeFromActiveBack( value );
   if( ptr == lastactive ) return removeFromBoundFront( value );
   if( tList< NodeType >::removeNext( value, ptr ) ){
      int act = value.getBoundaryFlag();
      if ((act == kNonBoundary) || (act == kStream )) nActiveNodes--;
      return 1;
   }
   return 0;
}

//SMM - added 08132008
template< class NodeType >
int tMeshList< NodeType >::
nextToBack( tListNode< NodeType > * prev )
{
   // Check if there was a previous node given and 
   // that there is a next node
   if( prev == 0 ) return 0;
   if( prev->next == 0 ) return 0;
   if (prev == this->lastactive) return 0;
   // This is the node that will be moved
   tListNode< NodeType > *nnode = prev->next;
   nActiveNodes--;
   // Hook up previous to the next node
   prev->next = nnode->next;
   // Check if last active and make previous so
   if (this->lastactive == nnode) this->lastactive = prev;
   // Hook up moved node to end
   this->last->next = nnode;
   nnode->next = 0;
   this->last = nnode;
   return 1;
}

//SMM - added 08132008
template< class NodeType >
int tMeshList< NodeType >::
frontToBack()
{
   // Check if list is empty
   if ( this->isEmpty() ) return 0;
   // This is the node to be moved
   tListNode< NodeType > *nnode = this->first;
   nActiveNodes--;
   // If lastactive, set to 0
   if (nnode == this->lastactive) this->lastactive = 0;
   // If no other node, done
   if (nnode->next == 0) return 1;
   // Make the next node first
   this->first = nnode->next;
   // Move node to end
   this->last->next = nnode;
   nnode->next = 0;
   this->last = nnode;
   return 1;
}

template< class NodeType >
int tMeshList< NodeType >::
removeFromBoundFront( NodeType &value )
{
   assert( &value != 0 );
   if( this->isEmpty() ) return 0;
   else if( this->last == lastactive ) return 0;
   else{
      tListNode< NodeType > * temp = lastactive->next;
      if( this->first == this->last ) this->first = this->last = 0;
      else lastactive->next = lastactive->next->next;
      value = temp->data;
      delete temp;
      this->nNodes--;
      return 1;
   }
}

template< class NodeType >
int tMeshList< NodeType >::
removeFromActiveBack( NodeType &value )
{
   if( this->isEmpty() ) return 0;
   else{
      tListNode< NodeType > * temp = lastactive;
      if( this->first == lastactive ) lastactive = 0;
      if( this->first == this->last ) this->first = this->last = 0;
      else{
         tListNode< NodeType > * current = this->first;
         while( current->next != lastactive ) current = current->next;
         current->next = lastactive->next;
         lastactive->next = 0;
         lastactive = current;
      }
      value = temp->data;
      delete temp;
      this->nNodes--;
      nActiveNodes--;
      return 1;
   }
}

#endif

//=========================================================================
//
//
//                      End of tMeshList.h
//
//
//=========================================================================
