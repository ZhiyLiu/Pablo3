#include <stdio.h>
#include <list>
#include <math.h>

#include "M3DFigure.h"
#include "InterfiguralConstraints.h"
#include "Geodesic.h"
#include "M3DPGAStats.h"
#include "M3DPGAPrimitiveStats.h"
#include "M3DInterpolater.h"

// dibyendu - cpns
#include "M3DCPNSStats.h"

#include "SimilarityTransform3D.h"
#include "M3DObject.h"
#include "M3DObjectFile.h"
#include "WorldSystem.h"
#include "Image3D.h"
#include <GL/gl.h>
#include <GL/glu.h>

#include <typeinfo>

//#define DEBUG
#define INVARIANTS


using namespace std;


M3DObject::M3DObject()
{
#ifdef BINARY
    name = new char[10];
    strcpy(name, "default");
#else
    name = NULL;
#endif
    nameCount = new int;
	*nameCount = 1;

	pga_stats = NULL;
	atom_pga_stats = NULL;

	// dibyendu - cpns
	cpns_stats = NULL ;

	type = NotAdaptive;
	transformation = NULL;
	that = NULL;
	wrld = NULL;

	flipped	= false;

	isDilated = false ;
	dilationFactorInModelUnits = 0.0 ;

	// This counter keeps track of variables pga_stats, atom_pga_stats, that, and wrld
	counter = new int;
	*counter = 1;

	////boundary = new SubdivBoundary();
}

M3DObject::M3DObject(const M3DObject & obj)
{
    M3DFigure * figurePtr,
              * newFigure;
    int numFigures,
        numFigureTrees,
        i;
    M3DFigureTreeNode * figureTree;

	// Copy the name
    if (obj.name != NULL)
    {
        name = new char[strlen(obj.name) + 1];
        strcpy(name, obj.name);
    }
    else
        name = NULL;
	nameCount = new int;
	*nameCount = 1;
	transformation = obj.transformation;

	// Note: NULL pointers are copied, when they exist.  While this presently
	// has no value, it may be necessary in the future to retain the "empty
	// slots".
    numFigures = obj.figures.size();
    for(i = 0; i < numFigures; i++)
    {
        figurePtr = obj.figures[i];
        if(figurePtr != NULL)
            newFigure = figurePtr->clone();		// Copies the figure name
        else
            newFigure = NULL;

        figures.insert(figures.end(), newFigure);
    }

    numFigureTrees = obj.getFigureTreeCount();
    for(i = 0; i < numFigureTrees; i++)
    {
        figureTree = obj.getFigureTreeRoot(i);
        if(figureTree != NULL)
            figureTrees.push_back(figureTree->copyPtr());
        else
            figureTrees.push_back(NULL);
    }

	////
	/*if (obj.getSubdivBoundary() != NULL)
		boundary = new SubdivBoundary(obj.getSubdivBoundary());
	else
		boundary = NULL;*/

	/*SubdivBoundary * newBoundary;
	for (i = 0; i < numFigures; i++)
	{
		if (obj.getSubdivBoundary(i) != NULL)
		{
			newBoundary = new SubdivBoundary(this, obj.getSubdivBoundary(i));
			boundaries.insert(boundaries.end(), newBoundary);
		}
		else
			newBoundary = NULL;
	}*/

	pga_stats = obj.pga_stats;
	atom_pga_stats = obj.atom_pga_stats;
	// dibyendu - cpns
	cpns_stats = obj.cpns_stats ;

	type = obj.type;
	that = obj.that;
	wrld = obj.wrld;
	flipped	= obj.flipped;

	isDilated = obj.isDilated ;
	dilationFactorInModelUnits = obj.dilationFactorInModelUnits ;

	counter = obj.counter;
    (*counter)++;
}

M3DObject & M3DObject::operator = (const M3DObject & obj)
{
    M3DFigure * figurePtr,
              * newFigure;
    int numFigures,
        numFigureTrees,
        i;
    M3DFigureTreeNode * figureTree;

	if (&obj == this)
		return *this;

    deleteAll();

	if (nameCount != obj.nameCount) {
		if( that && that != this ) {
			if( that->nameCount == nameCount && *nameCount == 2 ) {
				--*nameCount;
			}
		}
		if (--*nameCount == 0) {
			if (name != NULL)
				delete [] name;
			delete nameCount;
		}
	    name = obj.name;
		nameCount = obj.nameCount;
		(*nameCount)++;
	}

	transformation = obj.transformation;

    numFigures = obj.figures.size();
    for(i = 0; i < numFigures; i++)
    {
        figurePtr = obj.figures[i];
        if(figurePtr != NULL)
            newFigure = figurePtr->assign();	// Doesn't copy the figure invariants
        else
            newFigure = NULL;	// NULL pointers are retained in case needed in the future

        figures.insert(figures.end(), newFigure);
    }

    numFigureTrees = obj.getFigureTreeCount();
    for(i = 0; i < numFigureTrees; i++)
    {
        figureTree = obj.getFigureTreeRoot(i);
        if(figureTree != NULL)
            figureTrees.push_back(figureTree->copyPtr());
        else
            figureTrees.push_back(NULL);
    }

	/*if (boundary != NULL)
		delete boundary;
	boundary = new SubdivBoundary(obj.getSubdivBoundary());*/

	/*SubdivBoundary * boundaryPtr, newBoundary;
    for(i = 0; i < numFigures; i++)
    {
		boundaryPtr = obj.getSubdivBoundary(i);
		if (boundaryPtr != NULL)
			newBoundary = new SubdivBoundary(this, boundaryPtr);
		else
			newBoundary = NULL;

		boundaries.insert(boundaries.end(), newBoundary);
	}*/

	if (counter != obj.counter) {
		bool counterBasedDelete	= false;
		if( that && that != this ) {
			if( that->counter == counter && *counter == 2 ) {
				delete that;
				that	= NULL;
			}
		}
		if(--*counter == 0 || counterBasedDelete) {
			if( pga_stats ) {
				delete pga_stats;
			}
			if( atom_pga_stats ) {
				delete atom_pga_stats;
			}
			if( cpns_stats )
				delete cpns_stats ;

			if( that && that != this ) {
				delete that;
			}
			if( wrld ) {
				delete wrld;
			}
			delete counter;
		}
		pga_stats = obj.pga_stats;
		atom_pga_stats = obj.atom_pga_stats;
		cpns_stats = obj.cpns_stats ;
		type = obj.type;
		that = obj.that;
		wrld = obj.wrld;
		flipped	= obj.flipped;
		counter = obj.counter;
		(*counter)++;
	}

	isDilated = obj.isDilated ;
	dilationFactorInModelUnits = obj.dilationFactorInModelUnits ;

    return (*this);
}

M3DObject::~M3DObject()
{	

    deleteAll();

// invariants were a good idea for saving large memory, but not here.
// this is an intentional memory leak because fixing it is too hard.
//  -GST 6-nov-2006
#ifdef INVARIANTS
	// Responsibility of getting rid of that belongs to the last
	// instance left.
	bool counterBasedDelete = false;
	if( that && that != this ) {
		if( that->counter == counter && *counter == 2 ) {
			delete that;
			that	= NULL;
		}
	}

	if (--*counter == 0) {
		if( pga_stats ) {
			delete pga_stats;
		}
		if( atom_pga_stats ) {
			delete atom_pga_stats;
		}
		if( cpns_stats )
			delete cpns_stats ;
		if( wrld ) {
			delete wrld;
		}
		if( that && that != this ) {
			delete that;
		}
		delete counter;
	}

	if (--*nameCount == 0) {
		if (name != NULL)
			delete [] name;
		delete nameCount;
	}

	// transformation doesn't belong to this class
#endif
}

void M3DObject::restore()
{
	*this = *that;
}

void M3DObject::print(int markedFigIndex, int markedPrimIndex, bool dump,
					  ostream & out) const
{
    int numFigures,
        i;

    if(name != NULL)
        out << "Object name: \"" << name << "\"\n";
    else
        out << "Unnamed Object\n";

	cout << "    This object is ";
	if (orientation())
		cout << "right";
	else
		cout << "left";
	cout << "-handed.\n";

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        out << "Figure " << i << ':' << endl;
        if(figures[i] != NULL) {
//			if (i == markedFigIndex)                AGG: Why is this commented out?
//				figures[i]->print(dump, out, markedPrimIndex);
//			else
			figures[i]->print(dump, out);
		}
        else
            out << "NULL\n";

        out << '\n';
    }
	out << flush;
}

void M3DObject::print_tree_nodes(std::ostream & out) const
{
    for(int i = 0; i < figureTrees.size(); i++)
    {
        out << "    Tree #" << i;
        if(figureTrees[i] != NULL) {
           out << " (figure " << figureTrees[i]->getFigureId() << "):\n";
           figureTrees[i]->print(out);
		}
		else out << " - Null figure\n";
    }
	out << endl;
}

void M3DObject::printTree(std::ostream & out) const
{
	out << "Forest:\n";
	print_tree_nodes(out);
    out << flush;
}

char * M3DObject::copyName() const
{
    char * retStr;

    if(name == NULL)
        return name;

    retStr = new char[strlen(name) + 1];
    strcpy(retStr, name);

    return retStr;
}

void M3DObject::setName(const char * newName)
{
	if (--*nameCount == 0) {
		if (name != NULL)
			delete [] name;
	}
	else {
		nameCount = new int;
		*nameCount = 0;
	}

    if (newName != NULL)
    {
        name = new char[strlen(newName) + 1];
        strcpy(name, newName);
    }
    else
        name = NULL;

    (*nameCount)++;
}

void M3DObject::setPGAStats(M3DPGAStats * pgaStats, bool keep) {
	if (pga_stats == pgaStats)
		return;
	if (! keep && pga_stats != NULL)
		delete pga_stats;
	pga_stats = pgaStats;
}

void M3DObject::setCPNSStats( M3DCPNSStats * cpnsStats, bool keep ) {
	if( cpns_stats == cpnsStats )
		return ;
	if( keep == 0 && cpns_stats != NULL )
		delete cpns_stats ;
	cpns_stats = cpnsStats ;

}

void M3DObject::setAtomPGAStats(M3DPGAPrimitiveStats * pgaStats,
	bool keep)
{
	if (atom_pga_stats == pgaStats)
		return;
	if (! keep && atom_pga_stats != NULL)
		delete atom_pga_stats;
	atom_pga_stats = pgaStats;
}

int M3DObject::getLandmarkCount() const
{
    int lmCount,
        numFigures,
        i;

    lmCount = 0;
    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
        lmCount += figures[i]->getLandmarkCount();

    return lmCount;
}

int M3DObject::getMarkedLandmark(int & figureId)
{
	int numFigures, l;

    numFigures = figures.size();
    for (figureId = 0; figureId < numFigures; figureId++) {
        l = figures[figureId]->getMarkedLandmark();
		if (l >= 0)
			return l;
	}
	figureId = -1;
	return -1;
}

int M3DObject::getPrimitiveCount() const
{
    int primCount,
        numFigures,
        i;

    primCount = 0;
    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
        primCount += figures[i]->getPrimitiveCount();

    return primCount;
}

M3DFigure * M3DObject::getFigurePtr(int figIndex) const
{
    if(figIndex >= 0 && figIndex < figures.size())
        return figures[figIndex];

    return NULL;
}

void M3DObject::setFigurePtr(int figIndex, M3DFigure * fig_in)
{
	delete figures[figIndex];
	figures[figIndex] = fig_in;
}


M3DFigure * M3DObject::getFigurePtr(char * figName) const
{
    for (int figIndex = 0; figIndex < figures.size(); figIndex++) {
		if (0 == strcmp(figures[figIndex]->getName(), figName))
			return figures[figIndex];
	}
    return NULL;
}

int M3DObject::getFigureIndex(M3DFigure * figurePtr) const
{
    if (figurePtr == NULL)
		return -1;

	for (int i = 0; i < figures.size(); i++)
		if (figurePtr == figures[i])
			return i;

	return -1;
}

M3DPrimitive * M3DObject::getPrimitivePtr(int primIndex) const
{
    int primCount,
        numFigures,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        primCount = figures[i]->getPrimitiveCount();

        if(primIndex < primCount)
            return figures[i]->getPrimitivePtr(primIndex);

        primIndex -= primCount;
    }

    return NULL;
}

// Given a index of an atom in the object, return the index of the
// figure and change the atom index to that within the figure.
int M3DObject::getFigureAtomIndexes(int & primIndex) const
{
    int primCount,
        numFigures,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        primCount = figures[i]->getPrimitiveCount();

        if(primIndex < primCount)
            return i;

        primIndex -= primCount;
    }

    return -1;
}

M3DPrimitive * M3DObject::getPrimitivePtr(int figIndex, int primIndex) const
{
    if(figIndex < 0 || figIndex >= figures.size())
        return NULL;

    return figures[figIndex]->getPrimitivePtr(primIndex);
}

// Adds a figure and makes a tree to contain it

void M3DObject::addFigure(M3DFigure * figPtr)
{
    M3DFigureTreeNode * newRoot;

    // Insert figure into array of figures
    figures.insert(figures.end(), figPtr);

    // Create a new tree with this figure as the root
    newRoot = new M3DFigureTreeNode;
    newRoot->setFigureId(figures.size() - 1);
    newRoot->setParent(NULL);

    // Insert into forest of trees
    figureTrees.insert(figureTrees.end(), newRoot);

#ifdef DEBUG
	printTree();
#endif
}

bool M3DObject::testForPotentialLoop(int parentId, int childId)
{
	M3DFigure * parent;
	M3DFigure * child;
	M3DFigureTreeNode * parentTree;
	M3DFigureTreeNode * childTree;

#ifdef DEBUG
	cout << "Test for loop: parent = " << parentId << ", child = " << childId << endl;
#endif
	parent = getFigurePtr(parentId);
	child = getFigurePtr(childId);
	if (parent == NULL || child == NULL)
		return false;

	parentTree = getFigureTreeNode(parentId);
	if (parentTree->getChildCount() == 0)
		return false;

	childTree = getFigureTreeNode(childId);
	if (childTree->getParent() == NULL)
		return false;

	if (parentTree->findNode(childId) == NULL)
		return false;

	return true;
}

/*  This function is used by addTree(), below.  It copies the contents of tree
	to newTree and copies its associated figures from *oldObject to the current
	object.  Tree should be a node in *oldObject and newTree should be a node in
	*this.  The figure numbers in any constraints will not be changed and so will
	be incorrect after this operation.  However, the correspondence array will be
	updated, so they can be adjusted after all figures have been copied.
*/
void M3DObject::copyTreeFigures(M3DFigureTreeNode * tree, M3DObject * oldObject,
								M3DFigureTreeNode * newTree, int * correspondence)
{
	M3DFigure * parent;
	M3DFigure * newFigure;
	int numChildren;
	int parentId, newFigureId;
	int i;

	if (tree == NULL)
		return;

	// Add the parent figure
	parentId = tree->getFigureId();
	parent = oldObject->getFigurePtr(parentId);
	newFigure = parent->clone();
    figures.insert(figures.end(), newFigure);
	newFigureId = getFigureCount() - 1;
	correspondence[parentId] = newFigureId;

	// Make the id in the parent node refer to the new object
	newTree->setFigureId(newFigureId);

    // Copy the remaining figures
	numChildren = tree->getChildCount();
	for (i = 0; i < numChildren; i++) {
		M3DFigureTreeNode * subTree = tree->getChild(i);
		M3DFigureTreeNode * newSubTree = newTree->getChild(i);
		if (subTree != NULL)
			copyTreeFigures(subTree, oldObject, newSubTree, correspondence);
	}
}

/*  This function adds a copy of a tree along with its associated figures
    from *oldObject to the current object.  Tree must be a node in *oldObject.
	The number of trees in the current object will increase by 1.  Array
	correspondence is used to return the mapping of the figure IDs of the
	added figures to their new values in the combined object.  This is
	needed to correct the constraints' figure IDs, after all new trees from
	*oldObject are added.  The size of correspondence should equal the
	number of figures in *oldObject.  Any single call to this function will
	only fill the entries of correspondence of the figures added in that call.
	If a figure with figureId = x is added, then correspondence[x] will
	contain the new figureId of that figure after returning from the function.
*/
void M3DObject::addTree(M3DFigureTreeNode * tree, M3DObject * oldObject,
						int * correspondence)
{
    M3DFigureTreeNode * newTree;
	M3DFigure * newFigure;
	M3DFigure * root;
	int numChildren;
	int rootId, newFigureId;
	int i;

	if (tree == NULL)
		return;

	// Add the root figure
	rootId = tree->getFigureId();
	root = oldObject->getFigurePtr(rootId);
	newFigure = root->clone();

    figures.insert(figures.end(), newFigure);
	newFigureId = getFigureCount() - 1;
	correspondence[rootId] = newFigureId;

	newTree = tree->copyPtr();	// Copy the tree
	newTree->setFigureId(newFigureId);

    // Add the remaining figures
	numChildren = tree->getChildCount();
	for (i = 0; i < numChildren; i++) {
		M3DFigureTreeNode * subTree = tree->getChild(i);
		M3DFigureTreeNode * newSubTree = newTree->getChild(i);
		if (subTree != NULL)
			copyTreeFigures(subTree, oldObject, newSubTree, correspondence);
	}

	// Add the new tree
	figureTrees.insert(figureTrees.end(), newTree);

#ifdef DEBUG
	printTree();
#endif
}

M3DFigure * M3DObject::removeFigure(int index)
{
    M3DFigure * figurePtr;
    M3DFigureTreeNode * treeNode;
    M3DFigureTreeNode * treeRoot;
    M3DFigureTreeNode * parentNode;
    M3DFigureTreeNode * childNode;
    int numFigures,
        numTrees,
        childCount;
	int governorId, figureId;
    int i, j;

    numFigures = figures.size();
    if(index < 0 || index >= numFigures)
        return NULL;

    figurePtr = figures[index];

	// Remove constraints to and from this figure and adjust figure
	// numbers affected by removal of a figure
	InterfiguralConstraints & ifc = figurePtr->constraints();
	ifc.clear();

	for (governorId = 0; governorId < numFigures; governorId++) {
		bool deleteFigure = false;
		M3DFigure *governor = getFigurePtr(governorId);
		InterfiguralConstraints & ifc = governor->constraints();
		for (int i = 0; i < ifc.size(); i++) {
			figureId = ifc.figure(i);
			if (figureId == index) 
				deleteFigure = true;
			else if (figureId > index)	// Adjust higher index constraints for figure removal
				(void) ifc.changeFigureId(figureId, figureId - 1);
		}
		if (deleteFigure)
			(void) ifc.deleteFigure(index);
	}

    // Remove from the array of figures
    figures.erase(figures.begin() + index);

    // Remove node index from the tree, moving children to their own trees
    numTrees = figureTrees.size();
    for(i = 0; i < numTrees; i++)
    {
        treeRoot = figureTrees[i];
        if(treeRoot == NULL)
            continue;

        if(treeRoot->getFigureId() == index)
        {
            figureTrees.erase(figureTrees.begin() + i);
            --numTrees;
            treeNode = treeRoot;
        }
        else
            treeNode = treeRoot->findNode(index);

        if(treeNode != NULL)
        {
            // Move child subtrees to their own trees
            childCount = treeNode->getChildCount();
            for(j = 0; j < childCount; j++) {
                // Always pop the front element
                figureTrees.insert(figureTrees.end(), treeNode->popChild(0));
				// Remove links from new root nodes
				int n = figureTrees.size() - 1;
				int l = figureTrees[n]->getLinkCount();
				for (int k = l - 1; k >= 0; k--)
					figureTrees[n]->removeLink(k);
				// Also, clear the blend parameters in new root
				figureTrees[n]->setBlendAmount(0.0);
				figureTrees[n]->setBlendExtent(0.0);
				// Mark the new root as unattached
				figureTrees[n]->setAttachmentMode(M3DFigureTreeNode::UNATTACHED);
			}

            // Detach node from its parent and delete it
            parentNode = treeNode->getParent();
            if(parentNode != NULL)
            {
                childCount = parentNode->getChildCount();
                for(j = 0; j < childCount; j++)
                {
                    childNode = parentNode->getChild(j);
                    if(childNode != NULL && childNode == treeNode)
                    {
                        parentNode->popChild(j);
                        break;
                    }
                }
            }

            delete treeNode;
        }
    }

	// Update figureId numbers whenever greater than index
    numTrees = figureTrees.size();	// Must do this again
    for (i = 0; i < numTrees; i++)
    {
        treeRoot = figureTrees[i];
        if(treeRoot == NULL)
            continue;

		treeRoot->decrement(index);
    }

#ifdef DEBUG
	printTree();
#endif

    return figurePtr;
}

void M3DObject::deleteFigure(int index)
{
    M3DFigure * figurePtr;

    figurePtr = removeFigure(index);

    if(figurePtr != NULL)
        delete figurePtr;
}

void M3DObject::replaceFigure(int index, M3DFigure * figurePtr)
{
    M3DFigure * deleteFigure;

    if(index < 0 || index >= figures.size())
        return;

    deleteFigure = figures[index];
    figures[index] = figurePtr;
    delete deleteFigure;
}

M3DFigureTreeNode * M3DObject::getFigureTreeNode(int figureId)
{
    int numTrees;
    int i;

    M3DFigureTreeNode * treeNode;

    numTrees = figureTrees.size();

    for(i = 0; i < numTrees; i++)
    {
        treeNode = figureTrees[i];

        if(treeNode == NULL)
            continue;

        treeNode = treeNode->findNode(figureId);
        if(treeNode != NULL)
            return treeNode;
    }

    return NULL;
}

/*  Attach the tree (or partial tree) starting at treeNode to the
    specified figure.  Note that the tree containing the subtree to
	be attached may be the same as the tree to which it will be
	attached, resulting in moving of a subtree.  This function only
	changes the parent, children and siblings members of the nodes
	involved.  This function assumes that treeNode and targetNode
	(i.e., figureId) are both part of the same object (*this). 
*/
bool M3DObject::attachFigureTreeNode(M3DFigureTreeNode * treeNode, int figureId)
{
    int numTrees;
    int treeId, i;
    M3DFigureTreeNode * targetNode;
    M3DFigureTreeNode * parentNode;
	M3DFigureTreeNode * node;

	targetNode = getFigureTreeNode(figureId);
    if (targetNode == NULL) {
		// This should never happen
        cout << "Error: trying to attach an empty figure tree." << endl;
        return false;
    }
	if (treeNode->findNode(figureId) != NULL) {
		// This would create a loop
		cout << "Error: trying to attach a figure tree to itself." << endl;
        return false;
	}

	numTrees = figureTrees.size();
	if (treeNode->getParent() == NULL) {

		// Figure to be attached is a tree root

		// Find the node on the list of trees
		for (treeId = 0; treeId < numTrees; treeId++) {
			if (treeNode == figureTrees[treeId])
				break;
		}
		if (treeId > numTrees) {
			// This should never happen
			cout << "Error: trying to attach to an unknown figure tree." << endl;
	        return false;
		}
		for (i = treeId + 1; i < numTrees; i++)
			figureTrees[i - 1] = figureTrees[i];
		figureTrees.pop_back();
	}
	else {

		// Figure to be attached is a subfigure of a tree

		// Find the index of the tree containing the node to be attached
		if (numTrees == 1)
			treeId = 0;
		else {
			for (treeId = 0; treeId < numTrees; treeId++) {
				node = figureTrees[treeId]->findNode(targetNode->getFigureId());
				if (node != NULL)
					break;	// node == targetNode
			}
			if (treeId > numTrees) {
				// This should never happen
				cout << "Error: trying to attach to an unknown figure tree." << endl;
				return false;
			}
		}
		// Detach the node from the tree
		parentNode = treeNode->getParent();
		cout << "Detaching subfigure " << treeNode->getFigureId() 
			<< " from figure " << parentNode->getFigureId() << endl;
		treeNode->setParent(NULL);
		for (i = 0; i < parentNode->getChildCount(); i++) {
			node = parentNode->getChild(i);
			if (node == treeNode) {
				parentNode->popChild(i);
				break;
			}
		}
	}

	// Now attach the tree to the target node
	targetNode->addChild(treeNode);
    return true;
}

/*  Add a tree (or partial tree) starting at treeNode to the current
    object.  No figures are added and the figure Id's in the tree
	are not adjusted.  The figures must exist, however, because the
	hinge bit in some atoms will be set to correspond to the linking
	defined in the node being added.
*/
void M3DObject::addFigureTree(M3DFigureTreeNode * treeNode)
{
    if (treeNode == NULL)
    {
        cout << "Error: trying to add NULL figure tree." << endl;
        return;
    }

    // Insert as root of a new tree
    figureTrees.insert(figureTrees.end(), treeNode);
    markHingeAtoms(treeNode);
	int id = verifyConnectivity(treeNode);
	if (id >= 0)
		cout << "Warning: Figure " << id
			<< " may be miss-oriented relative to its main figure" << endl;
}

void M3DObject::markHingeAtoms(M3DFigureTreeNode * node)
{
    int i, j, count, nlinks;
    M3DFigure * figure;
	M3DPrimitive * atom;
    M3DFigureTreeNode * childNode;
	M3DPrimitiveLinkInfo * link;

	count = node->getChildCount();
	if (count == 0)
		return;

	for (i = 0; i < count; i++) {
		childNode = node->getChild(i);
		figure = getFigurePtr(childNode->getFigureId());
		nlinks = childNode->getLinkCount();
		for (j = 0; j < nlinks; j++) {
			// Set the hinge bit in each link atom
			link = childNode->getLink(j);
			atom = figure->getPrimitivePtr(link->primitiveId);
			atom->toggleHinge(true);
		}
		markHingeAtoms(childNode);
	}
}

// This function can be removed once it is certain that multifigure objects do
// not exist with incorrect linking.  If the model is verified to be correct, this
// function return -1.  Otherwise, it returns the figure Id of the first figure
// found to be misattached.
int M3DObject::verifyConnectivity(M3DFigureTreeNode * node)
{
    int i, j, count, nlinks;
	int figureId;
    M3DFigureTreeNode * childNode;
	M3DPrimitiveLinkInfo * link;
	int r, c, id;

	count = node->getChildCount();
	if (count == 0)
		return -1;

	for (i = 0; i < count; i++) {
		childNode = node->getChild(i);
		figureId = childNode->getFigureId();
		const M3DFigure* figure = getFigurePtr(figureId);
		if( typeid(*figure) == typeid(M3DTubeFigure) ) {
			r	= 1;
			c	= dynamic_cast<const M3DTubeFigure*>(figure)->getColumnCount();
		}
		else if( typeid(*figure) == typeid(M3DQuadFigure) ) {
			const M3DQuadFigure* fig	= dynamic_cast<const M3DQuadFigure*>(figure);
			r = fig->getRowCount();
			c = fig->getColumnCount();
		}
		else {
			// FIXME:
			cout << __FILE__ << ":" << __LINE__ << ": Unknown figure of type " << typeid(*figure).name() << endl;
			assert(false);
		}
		nlinks = childNode->getLinkCount();

		// Verify that the subfigure's row index increases away from the main figure
		for (j = 0; j < nlinks; j++) {
			link = childNode->getLink(j);
			id = link->primitiveId;
			if (id != c*j)
				return figureId;
		}
		return verifyConnectivity(childNode);
	}
	return -1;
}

/*  Detach a subtree and return a pointer to it.  If it cannot
    be detached, NULL is returned.  The blend settings are retained.
*/
M3DFigureTreeNode * M3DObject::detachFigureTreeNode(int figureId)
{
    int numTrees;
    int count;
    int i;
    M3DFigureTreeNode * rootNode;
    M3DFigureTreeNode * node;
    M3DFigureTreeNode * parentNode;
    M3DFigureTreeNode * childNode;
    M3DFigure * figure;
	M3DPrimitive * atom;

    numTrees = figureTrees.size();
	node = NULL;
	for (i = 0; i < numTrees; i++) {
            rootNode = figureTrees[i];
            if (rootNode == NULL)
                continue;
			node = rootNode->findNode(figureId);
            if (node != NULL)
				break;
	}
	if (node != NULL) {
		parentNode = node->getParent();
		if (parentNode == NULL)
			return NULL;	// Cannot detach a root node

		count = parentNode->getChildCount();
		for (i = 0; i < count; i++) {
			childNode = parentNode->getChild(i);
			if (childNode == node) {
				parentNode->popChild(i);
				break;
			}
		}

		figure = getFigurePtr(figureId);
		node->setParent(NULL);
		count = node->getLinkCount();
		for (i = count - 1; i >= 0; i--) {
			// Clear the hinge bit in each link atom
			M3DPrimitiveLinkInfo * link = node->getLink(i);
			atom = figure->getPrimitivePtr(link->primitiveId);
			atom->toggleHinge(false);

			node->removeLink(i);
		}
		figureTrees.push_back(node);	// Add back the detached tree

		node->setAttachmentMode(M3DFigureTreeNode::UNATTACHED);

		return node;
	}

	return NULL;
}


//
// Morphological operations:
// The units for the dilationFactor are in model co-ordinates
// i.e. 0-1 space
// dibyendu - modified this function to handle s-reps

void M3DObject::dilate( const double dilationFactor )
{	
	M3DPrimitive* prim;
	for( int i = 0; i != getPrimitiveCount(); ++i ) {
		prim	= getPrimitivePtr(i);		

		// separate dilation factor for quad primitives
		if( (dynamic_cast <M3DQuadPrimitive *> (prim)) != NULL )
			(dynamic_cast <M3DQuadPrimitive *> (prim))->dilate( dilationFactor ) ;

		prim->setR(prim->getR() + dilationFactor);
	}

	setDilationInfo( dilationFactor ) ;

}


void M3DObject::erode( const double erosionFactor )
{
	dilate(-erosionFactor);
}

void M3DObject::setDilationInfo(double _dilationFactorInModelUnits)
{
	// add the total dilation of the model
	dilationFactorInModelUnits += _dilationFactorInModelUnits ;

	// if total dilation factor is zero, then object is not dilated
	if( ( dilationFactorInModelUnits < 1e-9 ) && ( dilationFactorInModelUnits > -1e-9 ) )
		isDilated = 0 ;
	else
		isDilated = 1 ;

}

void M3DObject::subdivide()
{
	for( int i = 0; i != getFigureCount(); ++i ) {
		getFigurePtr(i)->subdivide();
	}
}


void M3DObject::resampleForRegularSpacing()
{
	for( int i = 0; i != getFigureCount(); ++i ) {
		getFigurePtr(i)->resampleForRegularSpacing();
	}
}


bool M3DObject::isRootNode(int figureId)
{
    M3DFigureTreeNode * node;

	node = getFigureTreeNode(figureId);
	if (node->getParent() == NULL)
		return true;
	else
		return false;
}

void M3DObject::select()
{
    vector<M3DFigure *>::iterator it;
    M3DFigure * figPtr;

    for(it = figures.begin(); it != figures.end(); it++)
    {
        figPtr = *it;
        if(figPtr != NULL)
            figPtr->select();
    }
}

void M3DObject::deselect()
{
    vector<M3DFigure *>::iterator it;
    M3DFigure * figPtr;

    for(it = figures.begin(); it != figures.end(); it++)
    {
        figPtr = *it;
        if(figPtr != NULL)
            figPtr->deselect();
    }
}

void M3DObject::toggleSelected()
{
    vector<M3DFigure *>::iterator it;
    M3DFigure * figPtr;

    for(it = figures.begin(); it != figures.end(); it++)
    {
        figPtr = *it;
        if(figPtr != NULL)
            figPtr->toggleSelected();
    }
}

bool M3DObject::isSelected()
{
    vector<M3DFigure *>::iterator it;
    M3DFigure * figPtr;

    for(it = figures.begin(); it != figures.end(); it++)
    {
        figPtr = *it;
        if(figPtr != NULL && (! figPtr->isSelected()))
            return false;
    }

    return true;
}

bool M3DObject::isAnySelected()
{
    vector<M3DFigure *>::iterator it;
    M3DFigure * figPtr;

    for(it = figures.begin(); it != figures.end(); it++)
    {
        figPtr = *it;
        if(figPtr != NULL && figPtr->isAnySelected())
            return true;
    }

    return false;
}

bool M3DObject::isModified()
{
    M3DFigure * figurePtr;
    int numFigures,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        figurePtr = figures[i];
        if(figurePtr != NULL && figurePtr->isModified())
            return true;
    }

    return false;
}

void M3DObject::setModified(bool flag)
{
    M3DFigure * figurePtr;
    int numFigures,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        figurePtr = figures[i];
        if(figurePtr != NULL)
            figurePtr->setModified(flag);
    }
}

Vector3D M3DObject::getCOG(bool all) const
{
    Vector3D center(0.0, 0.0, 0.0);
    M3DFigure * currFigurePtr;
    M3DPrimitive * currPrimitivePtr;
    int numFigures,
        numPrimitives,
        primCount,
        i, j;

    numFigures = figures.size();
    primCount = 0;

	if (all) {
		for(i = 0; i < numFigures; i++)
		{
			currFigurePtr = figures[i];

			if(currFigurePtr == NULL)
				continue;

			numPrimitives = currFigurePtr->getPrimitiveCount();
			for(j = 0; j < numPrimitives; j++)
			{
				currPrimitivePtr = currFigurePtr->getPrimitivePtr(j);
				if(currPrimitivePtr == NULL)
					continue;

				center += currPrimitivePtr->getX();
				primCount++;
			}
		}
	}
	else {
		for(i = 0; i < numFigures; i++)
		{
			currFigurePtr = figures[i];

			if(currFigurePtr == NULL)
				continue;

			numPrimitives = currFigurePtr->getPrimitiveCount();
			for(j = 0; j < numPrimitives; j++)
			{
				currPrimitivePtr = currFigurePtr->getPrimitivePtr(j);
				if(currPrimitivePtr == NULL)
					continue;

				if(currPrimitivePtr->isSelected())
				{
					center += currPrimitivePtr->getX();
					primCount++;
				}
			}
		}
	}

    if(primCount != 0)
        center /= (double) primCount;

    return center;
}

void M3DObject::scaleWidth(double scaleFactor)
{
    int numFigures,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        if(figures[i] != NULL)
            figures[i]->scaleWidth(scaleFactor);
    }
}

void M3DObject::scaleBy(double scaleFactor)
{
    scaleBy(scaleFactor, getCOG());
}

void M3DObject::scaleBy(double scaleFactor, const Vector3D & center)
{
    int numFigures,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        if(figures[i] != NULL)
            figures[i]->scaleBy(scaleFactor, center);
    }
}

void M3DObject::translateBy(const Vector3D & vTrans, bool selectAll)
{
    int numFigures,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        if(figures[i] != NULL)
            figures[i]->translateBy(vTrans, selectAll);
    }
}

void M3DObject::rotateBy(const Quat &q)
{
    rotateBy(q, getCOG());
}

void M3DObject::rotateBy(const Quat &q, const Vector3D & center)
{
    int numFigures,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        if(figures[i] != NULL)
            figures[i]->rotateBy(q, center);
    }
}

void M3DObject::applySimilarity(const SimilarityTransform3D & transform,
                                const Vector3D & center)
{
    int numFigures, i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        if(figures[i] != NULL)
            figures[i]->applySimilarity(transform, center);
    }
}

void M3DObject::applySimilarityAboutCOG(const SimilarityTransform3D * transform)
{
    applySimilarity(*transform, getCOG());
}

void M3DObject::deleteAll()
{
    int numFigures,
        numFigureTrees,
        i;

    numFigures = figures.size();
    for(i = 0; i < numFigures; i++)
    {
        if(figures[i] != NULL) {
            delete figures[i];
			figures[i] = NULL;
		}
    }

    figures.clear();

    numFigureTrees = figureTrees.size();
    for(i = 0; i < numFigureTrees; i++)
    {
        if(figureTrees[i] != NULL) {
            delete figureTrees[i];
			figureTrees[i] = NULL;
		}
    }

    figureTrees.clear();

	/*if (boundary != NULL)
		delete boundary;
	boundary = NULL;*/

	/*for(i = 0; i < numFigures; i++)
    {
		if (boundaries[i] != NULL)
		{
			delete boundaries[i];
			boundaries[i] = NULL;
		}
	}

	boundaries.clear();*/

}

// This function inverts the lists of constraints
void M3DObject::invertConstraints() {
	int num_figures;
	int governorId, figureId;

	num_figures = getFigureCount();
	for (governorId = 0; governorId < num_figures; governorId++) {
		M3DFigure *governor = getFigurePtr(governorId);
		InterfiguralConstraints & inverse = governor->inverseConstraints();
		inverse.clear();
	}

	for (governorId = 0; governorId < num_figures; governorId++) {
		M3DFigure *governor = getFigurePtr(governorId);
		InterfiguralConstraints & ifc = governor->constraints();
		for (int i = 0; i < ifc.size(); i++) {
			figureId = ifc.figure(i);
			double dist = ifc.distance(i);
			M3DFigure *figure = getFigurePtr(figureId);
			InterfiguralConstraints & inverse = figure->inverseConstraints();
			inverse.addFigure(governorId, dist);
		}
	}
}

// Make certain that all primitives are inside the unit cube
bool M3DObject::verifyInBounds() const {
	int num_figures;
	int figureId;

	num_figures = getFigureCount();
	for (figureId = 0; figureId < num_figures; figureId++) {
        if (figures[figureId] != NULL) {
            if (! figures[figureId]->verifyInBounds())
				return false;
		}
	}
	return true;
}

void M3DObject::renumberSubTree(M3DFigureTreeNode * node, int & newFigNum)
{
		node->setFigureId(newFigNum++);
		int nChildren = node->getChildCount();
		for (int c = 0; c < nChildren; c++) {
			M3DFigureTreeNode * child = node->getChild(c);
			if (child == NULL)
				continue;
			renumberSubTree(child, newFigNum);
		}
}

// Rearrange the storage of figures (change figure numbers) so they are
// in hierarchical, depth-first order.  The argument is the new list of
// figure numbers in depth-first order.
void M3DObject::renumber(int * newFigureNums)
{
    std::vector<M3DFigure *> newFigures;
    std::vector<M3DFigureTreeNode *> newFigureTrees;
	int newFigNum, oldFigNum;
	M3DFigure * figure;
	M3DFigureTreeNode * oldNode;
	M3DFigureTreeNode * newNode;

	for (newFigNum = 0; newFigNum < getFigureCount(); newFigNum++) {
		oldFigNum = newFigureNums[newFigNum];
		figure = getFigurePtr(oldFigNum);
		newFigures.push_back(figure->clone());

		oldNode = getFigureTreeNode(oldFigNum);
		if (oldNode && oldNode->getParent() == NULL) {
			// Renumber a whole tree at one time, starting from the root
			newNode = oldNode->copyPtr();
			newNode->setFigureId(newFigNum);
			int figureId = newFigNum + 1;
			int nChildren = newNode->getChildCount();
			for (int c = 0; c < nChildren; c++) {
				M3DFigureTreeNode * child = newNode->getChild(c);
				if (child == NULL)
					continue;
				renumberSubTree(child, figureId);
			}
			newFigureTrees.push_back(newNode);
		}
	}
	deleteAll();
	figures = newFigures;
	figureTrees = newFigureTrees;
}

#ifndef BINARY

// Regularize the model by adjusting the atom positions
void M3DObject::regularize(double stepsize, int iterations, bool verbose)
{
	cout << __FILE__ << ":" << __LINE__ << "M3DObject::regularize is deprecated, do not use." << endl;
#if 0
    int figureId;
    M3DQuadFigure * figurePtr;
    M3DQuadFigure * tmpFigure;
	Geodesic geo;

	// By default: all the atoms are treated as internal (non-crest) ones

	for (figureId = 0; figureId < getFigureCount(); figureId++)
	{
		if (verbose) cout << "Figure " << figureId << ": ";
		figurePtr = dynamic_cast<M3DQuadFigure*>( getFigurePtr(figureId));
		if (figurePtr == NULL)
			return;
		int	i,		
			r,
			c;
		M3DPrimitive	*primPtr,
						*m[4],
						mTmp,
						mTmp2;
		Quat qTmp;
		Vector3D		vTmp,
						bTmp;

		for (i = 0; i < figurePtr->getPrimitiveCount(); i++)
		{
			primPtr = figurePtr->getPrimitivePtr(i);
			if (primPtr->getQ().getW() < 0)
			{
				qTmp = primPtr->getQ();
				qTmp.neg();
				primPtr->setQ(qTmp);
			}
		}

		tmpFigure = new M3DQuadFigure(*figurePtr);
		int rNum = tmpFigure->getRowCount(),
			cNum = tmpFigure->getColumnCount();

		for (i = 0; i < iterations; i++)
		{
			if (verbose) cout << ".";
			for (r = 1; r < rNum - 1; r++)
				for (c = 1; c < cNum - 1; c++)
				{
					primPtr = tmpFigure->getPrimitivePtr(r, c);
					m[0] = figurePtr->getPrimitivePtr(r-1, c);
					m[1] = figurePtr->getPrimitivePtr(r, c-1);
					m[2] = figurePtr->getPrimitivePtr(r+1, c);
					m[3] = figurePtr->getPrimitivePtr(r, c+1);

					geo.atomAverage(4, m, &mTmp);
#ifdef Debug
					cout << "average Atom for r=" << r << ", c=" << c <<endl;
#endif

					geo.atomInterp(stepsize, primPtr, &mTmp, &mTmp2);

					mTmp2.setR(primPtr->getR());

					vTmp = mTmp2.getX()-primPtr->getX();
					vTmp = vTmp-primPtr->getN()*(vTmp*primPtr->getN());

					mTmp2.setX(primPtr->getX()+vTmp);

					*primPtr = mTmp2;
				}

			*figurePtr = *tmpFigure;
		}
		if (verbose) cout << endl;

		delete tmpFigure;
	}
#endif
}

#else	/* BINARY */

double approximateT(Vector3D v, Vector3D deriv)
{
	Vector3D	t;
	int			devid=3;

	if(deriv.getX()!=0)
		t.setX(v.getX()/deriv.getX());
	else
	{
		t.setX(0);
		devid--;
	}
	if(deriv.getY()!=0)
		t.setY(v.getY()/deriv.getY());
	else
	{
		t.setY(0);
		devid--;
	}
	if(deriv.getZ()!=0)
		t.setZ(v.getZ()/deriv.getZ());
	else
	{
		t.setZ(0);
		devid--;
	}
	if(devid>0)
		return (t.getX()+t.getY()+t.getZ())/devid;
	else
		return 0;
}

// Regularize the model by adjusting the atom positions
void M3DObject::regularize(double stepsize, int iterations, bool verbose)
{
	// FIXME
	cerr << "THIS FUNCTION NEEDS A WORK-OVER due to tubes - rohit. Don't use it\n";
	/*
	 * This function uses the old Geodesic class, I wonder if this does anything
	 * useful right now - rrs.
	 * I haven't seen it being used in any case.
	 */
#if 0
	int figureId;
    M3DQuadFigure	*figurePtr,
					*tmpFigure;

	// by default: all the atoms are treated as standard (internal) ones
	Geodesic geo;
	// to deal with boundary atoms
	Geodesic geoE(false);

	for(figureId = 0; figureId<getFigureCount(); figureId++)
	{
		if (verbose) cout << "Figure " << figureId << ": ";
	figurePtr = dynamic_cast<M3DQuadFigure *>(getFigurePtr(figureId));
		if(figurePtr == NULL)
			return;
		int	i,		
			r,
			c;

		// for atom regularization
		M3DPrimitive	*primPtr;

		// for internal atoms
	M3DPrimitive	*m[4];
	M3DQuadPrimitive mTmp,
						mTmp2,
						p[3][3], pp[3][3];

		// for boundary atoms
	M3DQuadEndPrimitive	mE[3],
						mTmpE,
						mTmp2E;
		double			tt;

		Quat qTmp;
		Vector3D		vTmp,
						bTmp,
						tTmp, tTmp2;

		tmpFigure = new M3DQuadFigure(*figurePtr);
		int rNum = tmpFigure->getRowCount(),
			cNum = tmpFigure->getColumnCount();

		for(i=0; i<iterations; i++)
		{
			if (verbose == 1) cout << ".";

			// boundary atoms
			//*
			for(r=1; r<rNum-1; r++)
			{
				c=0;
				primPtr=tmpFigure->getPrimitivePtr(r, c);
				m[0]=figurePtr->getPrimitivePtr(r-1, c);
				m[1]=figurePtr->getPrimitivePtr(r+1, c);
				//m[2]=figurePtr->getPrimitivePtr(r, c+1);
				geoE.atomAverage(2, m, &mTmpE);
//#define debug_reg
#ifdef debug
				cout << "average Atom for r=" << r << ", c=" << c << " :" << endl;
				mTmpE.print();
#endif
				geoE.atomInterp(stepsize, primPtr, &mTmpE, &mTmp2E);
			//mE[1]=*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr)); working right here
			mE[0]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[0]));
			mE[2]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[1]));
				mE[1].setX(2*primPtr->getX()-0.5*mE[0].getX()-0.5*mE[2].getX());
				mE[1].setR(2*primPtr->getR()-0.5*mE[0].getR()-0.5*mE[2].getR());
			mE[1].setElongation(2*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation()-
									0.5*mE[0].getElongation()-0.5*mE[2].getElongation());
				tTmp=mE[2].getX()-mE[0].getX();
				tTmp2=tTmp;
				tTmp.normalize();
				// first to project deltaX to curve tangent
				vTmp=mTmp2E.getX()-primPtr->getX();
				vTmp=tTmp*(tTmp*vTmp);
				tt=.5+approximateT(vTmp, tTmp2);
				mTmp2E.setR((1-tt)*(1-tt)*mE[0].getR()+2*tt*(1-tt)*mE[1].getR()+tt*tt*mE[2].getR());
				mTmp2E.setElongation((1-tt)*(1-tt)*mE[0].getElongation()+2*tt*(1-tt)*mE[1].getElongation()+tt*tt*mE[2].getElongation());
				// then to the tangent plane
				vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
				mTmp2E.setX(vTmp+primPtr->getX());
				// 01/08/2004
				//mTmp2E.setR(primPtr->getR());
			//mTmp2E.setElongation((dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation());
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->setElongation(mTmp2E.getElongation());
				primPtr->setX(mTmp2E.getX());
				primPtr->setR(mTmp2E.getR());
				primPtr->setQ(mTmp2E.getQ());
				primPtr->setTheta(mTmp2E.getTheta());
#ifdef debug_reg
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->print();
#endif

				c=cNum-1;
				primPtr=tmpFigure->getPrimitivePtr(r, c);
				m[0]=figurePtr->getPrimitivePtr(r-1, c);
				m[1]=figurePtr->getPrimitivePtr(r+1, c);
				//m[2]=figurePtr->getPrimitivePtr(r, c-1);
				geoE.atomAverage(2, m, &mTmpE);
#ifdef debug
				cout << "average Atom for r=" << r << ", c=" << c << " :" << endl;
				mTmpE.print();
#endif
				geoE.atomInterp(stepsize, primPtr, &mTmpE, &mTmp2E);
			//mE[1]=*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr)); working right here
			mE[0]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[0]));
			mE[2]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[1]));
				mE[1].setX(2*primPtr->getX()-0.5*mE[0].getX()-0.5*mE[2].getX());
				mE[1].setR(2*primPtr->getR()-0.5*mE[0].getR()-0.5*mE[2].getR());
			mE[1].setElongation(2*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation()-
									0.5*mE[0].getElongation()-0.5*mE[2].getElongation());
				tTmp=mE[2].getX()-mE[0].getX();
				tTmp2=tTmp;
				tTmp.normalize();
				// first to project deltaX to curve tangent
				vTmp=mTmp2E.getX()-primPtr->getX();
				vTmp=tTmp*(tTmp*vTmp);
				tt=.5+approximateT(vTmp, tTmp2);
				mTmp2E.setR((1-tt)*(1-tt)*mE[0].getR()+2*tt*(1-tt)*mE[1].getR()+tt*tt*mE[2].getR());
				mTmp2E.setElongation((1-tt)*(1-tt)*mE[0].getElongation()+2*tt*(1-tt)*mE[1].getElongation()+tt*tt*mE[2].getElongation());
				// then to the tangent plane
				vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
				mTmp2E.setX(vTmp+primPtr->getX());
				// 01/08/2004
				//mTmp2E.setR(primPtr->getR());
			//mTmp2E.setElongation((dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation());
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->setElongation(mTmp2E.getElongation());
				primPtr->setX(mTmp2E.getX());
				primPtr->setR(mTmp2E.getR());
				primPtr->setQ(mTmp2E.getQ());
				primPtr->setTheta(mTmp2E.getTheta());
#ifdef debug_reg
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->print();
#endif
			}
			for(c=1; c<cNum-1; c++)
			{
				r=0;
				primPtr=tmpFigure->getPrimitivePtr(r, c);
				m[0]=figurePtr->getPrimitivePtr(r, c-1);
				m[1]=figurePtr->getPrimitivePtr(r, c+1);
				//m[2]=figurePtr->getPrimitivePtr(r+1, c);
				geoE.atomAverage(2, m, &mTmpE);
#ifdef debug
				cout << "average Atom for r=" << r << ", c=" << c <<endl;
#endif
				geoE.atomInterp(stepsize, primPtr, &mTmpE, &mTmp2E);
			//mE[1]=*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr)); working right here
			mE[0]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[0]));
			mE[2]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[1]));
				mE[1].setX(2*primPtr->getX()-0.5*mE[0].getX()-0.5*mE[2].getX());
				mE[1].setR(2*primPtr->getR()-0.5*mE[0].getR()-0.5*mE[2].getR());
			mE[1].setElongation(2*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation()-
									0.5*mE[0].getElongation()-0.5*mE[2].getElongation());
				tTmp=mE[2].getX()-mE[0].getX();
				tTmp2=tTmp;
				tTmp.normalize();
				// first to project deltaX to curve tangent
				vTmp=mTmp2E.getX()-primPtr->getX();
				vTmp=tTmp*(tTmp*vTmp);
				tt=.5+approximateT(vTmp, tTmp2);
				mTmp2E.setR((1-tt)*(1-tt)*mE[0].getR()+2*tt*(1-tt)*mE[1].getR()+tt*tt*mE[2].getR());
				mTmp2E.setElongation((1-tt)*(1-tt)*mE[0].getElongation()+2*tt*(1-tt)*mE[1].getElongation()+tt*tt*mE[2].getElongation());
				// then to the tangent plane
				vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
				mTmp2E.setX(vTmp+primPtr->getX());
				// 01/08/2004
				//mTmp2E.setR(primPtr->getR());
			//mTmp2E.setElongation((dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation());
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->setElongation(mTmp2E.getElongation());
				primPtr->setX(mTmp2E.getX());
				primPtr->setR(mTmp2E.getR());
				primPtr->setQ(mTmp2E.getQ());
				primPtr->setTheta(mTmp2E.getTheta());
#ifdef debug_reg
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->print();
#endif

				r=rNum-1;
				primPtr=tmpFigure->getPrimitivePtr(r, c);
				m[0]=figurePtr->getPrimitivePtr(r, c-1);
				m[1]=figurePtr->getPrimitivePtr(r, c+1);
				//m[2]=figurePtr->getPrimitivePtr(r-1, c);
				geoE.atomAverage(2, m, &mTmpE);
#ifdef debug
				cout << "average Atom for r=" << r << ", c=" << c <<endl;
#endif
				geoE.atomInterp(stepsize, primPtr, &mTmpE, &mTmp2E);
			//mE[1]=*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr)); working right here
			mE[0]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[0]));
			mE[2]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[1]));
				mE[1].setX(2*primPtr->getX()-0.5*mE[0].getX()-0.5*mE[2].getX());
				mE[1].setR(2*primPtr->getR()-0.5*mE[0].getR()-0.5*mE[2].getR());
			mE[1].setElongation(2*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation()-
									0.5*mE[0].getElongation()-0.5*mE[2].getElongation());
				tTmp=mE[2].getX()-mE[0].getX();
				tTmp2=tTmp;
				tTmp.normalize();
				// first to project deltaX to curve tangent
				vTmp=mTmp2E.getX()-primPtr->getX();
				vTmp=tTmp*(tTmp*vTmp);
				tt=.5+approximateT(vTmp, tTmp2);
				mTmp2E.setR((1-tt)*(1-tt)*mE[0].getR()+2*tt*(1-tt)*mE[1].getR()+tt*tt*mE[2].getR());
				mTmp2E.setElongation((1-tt)*(1-tt)*mE[0].getElongation()+2*tt*(1-tt)*mE[1].getElongation()+tt*tt*mE[2].getElongation());
				// then to the tangent plane
				vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
				mTmp2E.setX(vTmp+primPtr->getX());
				// 01/08/2004
				//mTmp2E.setR(primPtr->getR());
			//mTmp2E.setElongation((dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation());
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->setElongation(mTmp2E.getElongation());
				primPtr->setX(mTmp2E.getX());
				primPtr->setR(mTmp2E.getR());
				primPtr->setQ(mTmp2E.getQ());
				primPtr->setTheta(mTmp2E.getTheta());
#ifdef debug_reg
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->print();
#endif
			}
			//*/

			*figurePtr=*tmpFigure;
		}

		// the corner atoms
		if(iterations>0)
		{
				primPtr=tmpFigure->getPrimitivePtr(0, 0);
				m[0]=figurePtr->getPrimitivePtr(1, 0);
				m[1]=figurePtr->getPrimitivePtr(0, 1);
				//m[2]=figurePtr->getPrimitivePtr(r-1, c);
				geoE.atomAverage(2, m, &mTmpE);
#ifdef debug
				cout << "average Atom for r=" << r << ", c=" << c <<endl;
#endif
				geoE.atomInterp(stepsize, primPtr, &mTmpE, &mTmp2E);
			//mE[1]=*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr)); working right here
			mE[0]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[0]));
			mE[2]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[1]));
				mE[1].setX(2*primPtr->getX()-0.5*mE[0].getX()-0.5*mE[2].getX());
				mE[1].setR(2*primPtr->getR()-0.5*mE[0].getR()-0.5*mE[2].getR());
			mE[1].setElongation(2*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation()-
									0.5*mE[0].getElongation()-0.5*mE[2].getElongation());
				tTmp=mE[2].getX()-mE[0].getX();
				tTmp2=tTmp;
				tTmp.normalize();
				// first to project deltaX to curve tangent
				vTmp=mTmp2E.getX()-primPtr->getX();
				vTmp=tTmp*(tTmp*vTmp);
				tt=.5+approximateT(vTmp, tTmp2);
				mTmp2E.setR((1-tt)*(1-tt)*mE[0].getR()+2*tt*(1-tt)*mE[1].getR()+tt*tt*mE[2].getR());
				mTmp2E.setElongation((1-tt)*(1-tt)*mE[0].getElongation()+2*tt*(1-tt)*mE[1].getElongation()+tt*tt*mE[2].getElongation());
				// then to the tangent plane
				vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
				mTmp2E.setX(vTmp+primPtr->getX());
				// 01/08/2004
				//mTmp2E.setR(primPtr->getR());
			//mTmp2E.setElongation((dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation());
//			*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))=mTmp2E;
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->setElongation(mTmp2E.getElongation());
				primPtr->setX(mTmp2E.getX());
				primPtr->setR(mTmp2E.getR());
				primPtr->setQ(mTmp2E.getQ());
				primPtr->setTheta(mTmp2E.getTheta());

				primPtr=tmpFigure->getPrimitivePtr(rNum-1, cNum-1);
				m[0]=figurePtr->getPrimitivePtr(rNum-2, cNum-1);
				m[1]=figurePtr->getPrimitivePtr(rNum-1, cNum-2);
				//m[2]=figurePtr->getPrimitivePtr(r-1, c);
				geoE.atomAverage(2, m, &mTmpE);
#ifdef debug
				cout << "average Atom for r=" << r << ", c=" << c <<endl;
#endif
				geoE.atomInterp(stepsize, primPtr, &mTmpE, &mTmp2E);
			//mE[1]=*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr)); working right here
			mE[0]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[0]));
			mE[2]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[1]));
				mE[1].setX(2*primPtr->getX()-0.5*mE[0].getX()-0.5*mE[2].getX());
				mE[1].setR(2*primPtr->getR()-0.5*mE[0].getR()-0.5*mE[2].getR());
			mE[1].setElongation(2*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation()-
									0.5*mE[0].getElongation()-0.5*mE[2].getElongation());
				tTmp=mE[2].getX()-mE[0].getX();
				tTmp2=tTmp;
				tTmp.normalize();
				// first to project deltaX to curve tangent
				vTmp=mTmp2E.getX()-primPtr->getX();
				vTmp=tTmp*(tTmp*vTmp);
				tt=.5+approximateT(vTmp, tTmp2);
				mTmp2E.setR((1-tt)*(1-tt)*mE[0].getR()+2*tt*(1-tt)*mE[1].getR()+tt*tt*mE[2].getR());
				mTmp2E.setElongation((1-tt)*(1-tt)*mE[0].getElongation()+2*tt*(1-tt)*mE[1].getElongation()+tt*tt*mE[2].getElongation());
				// then to the tangent plane
				vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
				mTmp2E.setX(vTmp+primPtr->getX());
				// 01/08/2004
				//mTmp2E.setR(primPtr->getR());
			//mTmp2E.setElongation((dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation());
//			*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))=mTmp2E;
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->setElongation(mTmp2E.getElongation());
				primPtr->setX(mTmp2E.getX());
				primPtr->setR(mTmp2E.getR());
				primPtr->setQ(mTmp2E.getQ());
				primPtr->setTheta(mTmp2E.getTheta());

				primPtr=tmpFigure->getPrimitivePtr(rNum-1, 0);
				m[0]=figurePtr->getPrimitivePtr(rNum-2, 0);
				m[1]=figurePtr->getPrimitivePtr(rNum-1, 1);
				//m[2]=figurePtr->getPrimitivePtr(r-1, c);
				geoE.atomAverage(2, m, &mTmpE);
#ifdef debug
				cout << "average Atom for r=" << r << ", c=" << c <<endl;
#endif
				geoE.atomInterp(stepsize, primPtr, &mTmpE, &mTmp2E);
			//mE[1]=*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr)); working right here
			mE[0]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[0]));
			mE[2]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[1]));
				mE[1].setX(2*primPtr->getX()-0.5*mE[0].getX()-0.5*mE[2].getX());
				mE[1].setR(2*primPtr->getR()-0.5*mE[0].getR()-0.5*mE[2].getR());
			mE[1].setElongation(2*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation()-
									0.5*mE[0].getElongation()-0.5*mE[2].getElongation());
				tTmp=mE[2].getX()-mE[0].getX();
				tTmp2=tTmp;
				tTmp.normalize();
				// first to project deltaX to curve tangent
				vTmp=mTmp2E.getX()-primPtr->getX();
				vTmp=tTmp*(tTmp*vTmp);
				tt=.5+approximateT(vTmp, tTmp2);
				mTmp2E.setR((1-tt)*(1-tt)*mE[0].getR()+2*tt*(1-tt)*mE[1].getR()+tt*tt*mE[2].getR());
				mTmp2E.setElongation((1-tt)*(1-tt)*mE[0].getElongation()+2*tt*(1-tt)*mE[1].getElongation()+tt*tt*mE[2].getElongation());
				// then to the tangent plane
				vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
				mTmp2E.setX(vTmp+primPtr->getX());
				// 01/08/2004
				//mTmp2E.setR(primPtr->getR());
			//mTmp2E.setElongation((dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation());
//			*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))=mTmp2E;
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->setElongation(mTmp2E.getElongation());
				primPtr->setX(mTmp2E.getX());
				primPtr->setR(mTmp2E.getR());
				primPtr->setQ(mTmp2E.getQ());
				primPtr->setTheta(mTmp2E.getTheta());

				primPtr=tmpFigure->getPrimitivePtr(0, cNum-1);
				m[0]=figurePtr->getPrimitivePtr(1, cNum-1);
				m[1]=figurePtr->getPrimitivePtr(0, cNum-2);
				//m[2]=figurePtr->getPrimitivePtr(r-1, c);
				geoE.atomAverage(2, m, &mTmpE);
#ifdef debug
				cout << "average Atom for r=" << r << ", c=" << c <<endl;
#endif
				geoE.atomInterp(stepsize, primPtr, &mTmpE, &mTmp2E);
			//mE[1]=*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr)); working right here
			mE[0]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[0]));
			mE[2]=*(dynamic_cast<M3DQuadEndPrimitive *>(m[1]));
				mE[1].setX(2*primPtr->getX()-0.5*mE[0].getX()-0.5*mE[2].getX());
				mE[1].setR(2*primPtr->getR()-0.5*mE[0].getR()-0.5*mE[2].getR());
			mE[1].setElongation(2*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation()-
									0.5*mE[0].getElongation()-0.5*mE[2].getElongation());
				tTmp=mE[2].getX()-mE[0].getX();
				tTmp2=tTmp;
				tTmp.normalize();
				// first to project deltaX to curve tangent
				vTmp=mTmp2E.getX()-primPtr->getX();
				vTmp=tTmp*(tTmp*vTmp);
				tt=.5+approximateT(vTmp, tTmp2);
				mTmp2E.setR((1-tt)*(1-tt)*mE[0].getR()+2*tt*(1-tt)*mE[1].getR()+tt*tt*mE[2].getR());
				mTmp2E.setElongation((1-tt)*(1-tt)*mE[0].getElongation()+2*tt*(1-tt)*mE[1].getElongation()+tt*tt*mE[2].getElongation());
				// then to the tangent plane
				vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
				mTmp2E.setX(vTmp+primPtr->getX());
				// 01/08/2004
				//mTmp2E.setR(primPtr->getR());
			//mTmp2E.setElongation((dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->getElongation());
//			*(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))=mTmp2E;		
			(dynamic_cast<M3DQuadEndPrimitive *>(primPtr))->setElongation(mTmp2E.getElongation());
				primPtr->setX(mTmp2E.getX());
				primPtr->setR(mTmp2E.getR());
				primPtr->setQ(mTmp2E.getQ());
				primPtr->setTheta(mTmp2E.getTheta());

		*figurePtr=*tmpFigure;
		}

	//*
		// internal atoms
		static int gridNum=7;
		Vector3D *gridPnt=new Vector3D[gridNum+3];

		for(i=0; i<iterations; i++)
		{
			// cout << ".";
			// internal atoms
			//*
			for(r=1; r<rNum-1; r++)
				for(c=1; c<cNum-1; c++)
				{
					primPtr=tmpFigure->getPrimitivePtr(r, c);
					m[0]=figurePtr->getPrimitivePtr(r-1, c);
					m[1]=figurePtr->getPrimitivePtr(r, c-1);
					m[2]=figurePtr->getPrimitivePtr(r+1, c);
					m[3]=figurePtr->getPrimitivePtr(r, c+1);

					geo.atomAverage(4, m, &mTmp);
#ifdef debug
					cout << "average Atom for r=" << r << ", c=" << c <<endl;
#endif

#define BSPLINE_AT_ALL
//#define ALTERNATIVE_BSPLINE
#ifdef BSPLINE_AT_ALL
					// the X=(x, y, z, r) will be re-interpolated by a local b-spline patch
				p[0][0]=*(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r-1, c-1)));
				pp[0][1]=*(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r-1, c)));
				p[0][2]=*(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r-1, c+1)));
				pp[1][0]=*(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r, c-1)));
				//p[1][1]=(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r, c)));		// to be solved
				pp[1][2]=*(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r, c+1)));
				p[2][0]=*(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r+1, c-1)));
				pp[2][1]=*(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r+1, c)));
				p[2][2]=*(dynamic_cast<M3DQuadPrimitive*>(figurePtr->getPrimitivePtr(r+1, c+1)));

					p[0][1].setX(2*pp[0][1].getX()-0.5*p[0][0].getX()-0.5*p[0][2].getX());
					p[2][1].setX(2*pp[2][1].getX()-0.5*p[2][0].getX()-0.5*p[2][2].getX());
					p[1][0].setX(2*pp[1][0].getX()-0.5*p[0][0].getX()-0.5*p[2][0].getX());
					p[1][2].setX(2*pp[1][2].getX()-0.5*p[0][2].getX()-0.5*p[2][2].getX());

					// the calculation of p[1][1] is in each of the methods

					static double MAX_RATIO=1.6;
					double dis1, dis2;
					bool debugTag=false;
#ifdef ALTERNATIVE_BSPLINE
					p[1][1].setX((p[0][1].getX()+p[2][1].getX()+p[1][0].getX()+p[1][2].getX())/4);
					/*
					pp[1][1].setX((pp[0][1].getX()+pp[2][1].getX()+pp[1][0].getX()+pp[1][2].getX())/4);
					p[1][1].setX(2*pp[1][1].getX()-0.5*pp[1][0].getX()-0.5*pp[1][2].getX());
					p[1][1].setX(2*p[1][1].getX()-0.5*p[0][1].getX()-0.5*p[2][1].getX());
					//*/
					tTmp=0.25*p[0][1].getX()+0.5*p[1][1].getX()+0.25*p[2][1].getX();
					tTmp=0.25*pp[1][0].getX()+0.5*tTmp+0.25*pp[1][2].getX();
					dis1=tTmp.distSquare(pp[1][0].getX()); dis2=tTmp.distSquare(pp[1][2].getX());
					double scale;
					static double CONFIDENCE=1;
					if(dis1/dis2<=1/MAX_RATIO || dis1/dis2>=MAX_RATIO)
					{
						debugTag=true;
#ifdef debug
						printf("%f, ", dis1/dis2);
#endif
						scale=(1-1/(1+sqrt(dis1/dis2)*CONFIDENCE))*0.5;
					}
					else
					{
						scale=0.25;
					}
					p[1][1].setX(p[1][0].getX()*scale);
					p[1][1].setX(p[1][1].getX()+p[1][2].getX()*(0.5-scale));

					dis1=tTmp.distSquare(pp[0][1].getX()); dis2=tTmp.distSquare(pp[2][1].getX());
					if(dis1/dis2<=1/MAX_RATIO || dis1/dis2>=MAX_RATIO)
					{
						if(debugTag)
							cout << "! (" << r << ", " << c << ") !"; //printf("! (%d, %d) !", r, c);
						scale=(1-1/(1+sqrt(dis1/dis2)*CONFIDENCE))*0.5;
					}
					else
					{
						scale=0.25;
					}
					p[1][1].setX(p[1][1].getX()+p[0][1].getX()*scale);
					p[1][1].setX(p[1][1].getX()+p[2][1].getX()*(0.5-scale));

					tTmp=0.25*p[0][1].getX() +0.5*p[1][1].getX()+0.25*p[2][1].getX();
					tTmp=0.25*pp[1][0].getX()+0.5*tTmp          +0.25*pp[1][2].getX();

#ifdef debug
					if(debugTag)
					{
						printf("%f, ", dis1/dis2);
					}
					dis1=tTmp.distSquare(pp[1][0].getX()); dis2=tTmp.distSquare(pp[1][2].getX());
					if(debugTag)
					{
						printf("after: %f, ", dis1/dis2);
					}
					dis1=tTmp.distSquare(pp[0][1].getX()); dis2=tTmp.distSquare(pp[2][1].getX());
					if(debugTag)
					{
						printf("%f\n", dis1/dis2);
					}
#endif
					// the new bspline generator goes here : using an approximate 'medial' point of the patch
#else
					bool tTmpRecalculated=false;
					double param[3], segDis[50], totalLength=0, curLength=0;
					double paramStart=0, paramStep=1.0/(gridNum+1.0);
					p[1][1].setX((p[0][1].getX()+p[2][1].getX()+p[1][0].getX()+p[1][2].getX())/4);
					/*
					pp[1][1].setX((pp[0][1].getX()+pp[2][1].getX()+pp[1][0].getX()+pp[1][2].getX())/4);
					p[1][1].setX(2*pp[1][1].getX()-0.5*pp[1][0].getX()-0.5*pp[1][2].getX());
					p[1][1].setX(2*p[1][1].getX()-0.5*p[0][1].getX()-0.5*p[2][1].getX());
					//*/
					tTmp=0.25*p[0][1].getX()+0.5*p[1][1].getX()+0.25*p[2][1].getX();
					tTmp=0.25*pp[1][0].getX()+0.5*tTmp+0.25*pp[1][2].getX();
					dis1=tTmp.distSquare(pp[1][0].getX()); dis2=tTmp.distSquare(pp[1][2].getX());
					if(dis1/dis2<=1/MAX_RATIO || dis1/dis2>=MAX_RATIO)
					{
						debugTag=true;
#ifdef debug
						printf("%f, ", dis1/dis2);
#endif
	//					scale=(1-1/(1+sqrt(dis1/dis2)*CONFIDENCE))*0.5;

						tTmpRecalculated=true;
						int ii;
						for(ii=0; ii<gridNum+2; ii++)
						{
							param[0]=(1-paramStart)*(1-paramStart);
							param[1]=(paramStart)  *(1-paramStart)*2;
							param[2]=(paramStart)  *(paramStart);
							pp[1][1].setX(0.25*p[0][1].getX()+0.5*p[1][1].getX()+0.25*p[2][1].getX());
							gridPnt[ii]=param[0]*pp[1][0].getX()+param[1]*pp[1][1].getX()+param[2]*pp[1][2].getX();
							paramStart+=paramStep;
						}
						segDis[0]=0;
						for(ii=1; ii<gridNum+2; ii++)
						{
							segDis[ii]=sqrt(gridPnt[ii].distSquare(gridPnt[ii-1]));
							totalLength+=segDis[ii];
						}
						ii=1;
						while(curLength<0.5*totalLength)
						{
							curLength+=segDis[ii];
							ii++;
						}
						paramStart=(ii-0.5)*paramStep;
						param[0]=(1-paramStart)*(1-paramStart);
						param[1]=(paramStart)  *(1-paramStart)*2;
						param[2]=(paramStart)  *(paramStart);
						pp[1][1].setX(0.25*p[0][1].getX()+0.5*p[1][1].getX()+0.25*p[2][1].getX());
						gridPnt[gridNum+2]=param[0]*pp[1][0].getX()+param[1]*pp[1][1].getX()+param[2]*pp[1][2].getX();
						tTmp=gridPnt[gridNum+2];

						pp[1][1].setX(2*gridPnt[gridNum+2]-0.5*pp[1][0].getX()-0.5*pp[1][2].getX());
						p[1][1].setX(2*pp[1][1].getX()-0.5*p[1][0].getX()-0.5*p[1][2].getX());
					}
					else
					{
						dis1=tTmp.distSquare(pp[0][1].getX()); dis2=tTmp.distSquare(pp[2][1].getX());
						if(dis1/dis2<=1/MAX_RATIO || dis1/dis2>=MAX_RATIO)
						{
							debugTag=true;
#ifdef debug
							printf("%f, ", dis1/dis2);
#endif

						tTmpRecalculated=true;
						int ii;
						for(ii=0; ii<gridNum+2; ii++)
						{
							param[0]=(1-paramStart)*(1-paramStart);
							param[1]=(paramStart)  *(1-paramStart)*2;
							param[2]=(paramStart)  *(paramStart);
							pp[1][1].setX(0.25*p[1][0].getX()+0.5*p[1][1].getX()+0.25*p[1][2].getX());
							gridPnt[ii]=param[0]*pp[0][1].getX()+param[1]*pp[1][1].getX()+param[2]*pp[2][1].getX();
							paramStart+=paramStep;
						}
						segDis[0]=0;
						for(ii=1; ii<gridNum+2; ii++)
						{
							segDis[ii]=gridPnt[ii].distSquare(gridPnt[ii-1]);
							totalLength+=segDis[ii];
						}
						ii=1;
						while(curLength<0.5*totalLength)
						{
							curLength+=segDis[ii];
							ii++;
						}
						paramStart=(ii-0.5)*paramStep;
						param[0]=(1-paramStart)*(1-paramStart);
						param[1]=(paramStart)  *(1-paramStart)*2;
						param[2]=(paramStart)  *(paramStart);
						pp[1][1].setX(0.25*p[1][0].getX()+0.5*p[1][1].getX()+0.25*p[1][2].getX());
						gridPnt[gridNum+2]=param[0]*pp[0][1].getX()+param[1]*pp[1][1].getX()+param[2]*pp[2][1].getX();
						tTmp=gridPnt[gridNum+2];

						pp[1][1].setX(2*gridPnt[gridNum+2]-0.5*pp[0][1].getX()-0.5*pp[2][1].getX());
						p[1][1].setX(2*pp[1][1].getX()-0.5*p[1][0].getX()-0.5*p[1][2].getX());
						}
					}
#ifdef debug
					if(tTmpRecalculated)
					{
						dis1=tTmp.distSquare(pp[1][0].getX()); dis2=tTmp.distSquare(pp[1][2].getX());
						if(debugTag)
						{
							printf("after: %f, ", dis1/dis2);
						}
						dis1=tTmp.distSquare(pp[0][1].getX()); dis2=tTmp.distSquare(pp[2][1].getX());
						if(debugTag)
						{
							printf("%f\n", dis1/dis2);
						}
					}
#endif
#endif
					mTmp.setX(tTmp);
#endif
					geo.atomInterp(stepsize, primPtr, &mTmp, &mTmp2);
					mTmp2.setR(primPtr->getR());

					// project the displacement vector to the tangent plane at current atom
					vTmp=mTmp2.getX()-primPtr->getX();
					vTmp=vTmp-primPtr->getN()*(vTmp*primPtr->getN());
					mTmp2.setX(primPtr->getX()+vTmp);
	//				*primPtr=mTmp2;
				primPtr->setX(mTmp2.getX());
				primPtr->setR(mTmp2.getR());
				primPtr->setQ(mTmp2.getQ());
				primPtr->setTheta(mTmp2.getTheta());
				}

			*figurePtr=*tmpFigure;
		}
	//*/
		delete []gridPnt;

		if (verbose) cout << "Regularized" << endl;
		delete tmpFigure;
	}
	setModified(true);
#endif
}

#endif	/* BINARY */

/*	This function performs the following, composed transform
		M_1 -> W_1 -> W_2 -> M_2,
	where M refers to model coordinates and W refers to world
	coordinates, generally in cm.  The model-to-world and
	world-to-model transformations are specified in Image3D.
	The world-to-world transformation is merely a change of
	origin; no scaling is necessary, since all patient worlds
	are measured in cm in X, Y, and Z.
*/
bool M3DObject::applyWorld(Image3D * image)
{
	Vector3D newOrigin;
	const Vector3D * origOrigin;
	int numFigures, numAtoms, i, j;
	M3DFigure * figure;
	M3DPrimitive * atom;
	bool outside;

	if (wrld == NULL || image == NULL)
		return false;

	if (wrld->status() == false)
		return false;

    newOrigin = image->getWorldOrigin();
	origOrigin = wrld->getOrigin();

	numFigures = getFigureCount();
	outside = false;
	for (i = 0; i < numFigures; i++) {
		figure = getFigurePtr(i);
		if (figure == NULL)
			continue;
		numAtoms = figure->getPrimitiveCount();
		for (j = 0; j < numAtoms; j++) {
			atom = figure->getPrimitivePtr(j);
			if (atom == NULL)
				continue;

			Vector3D coords = atom->getX();
			wrld->modelToWorldCoordinates(coords);

			// World 1 to world 2 transformation
			coords -= *origOrigin;
			coords += newOrigin;

			image->worldToModelCoordinates(coords);
			atom->setX(coords);
			// Note: It is possible that the object lies off the
			// image and still is inside the unit cube.
			if (coords.getX() > 1.0 || coords.getX() < 0.0 ||
				coords.getY() > 1.0 || coords.getY() < 0.0 ||
				coords.getZ() > 1.0 || coords.getZ() < 0.0)
					outside = true;
		}
	}
	// Note that wrld is not changed to the new world coordinate space

	if (outside)
		cout << "Warning: model has been moved outside the unit cube" << endl;
	return true;
}

// The inverse of the previous function
bool M3DObject::unApplyWorld(Image3D * image)
{
	const Vector3D * newOrigin;
	Vector3D origOrigin;
	int numFigures, numAtoms, i, j;
	M3DFigure * figure;
	M3DPrimitive * atom;

	if (wrld == NULL || image == NULL)
		return false;

	if (wrld->status() == false)
		return false;

	newOrigin = wrld->getOrigin();
    origOrigin = image->getWorldOrigin();

	numFigures = getFigureCount();
	for (i = 0; i < numFigures; i++) {
		figure = getFigurePtr(i);
		if (figure == NULL)
			continue;
		numAtoms = figure->getPrimitiveCount();
		for (j = 0; j < numAtoms; j++) {
			atom = figure->getPrimitivePtr(j);
			if (atom == NULL)
				continue;

			Vector3D coords = atom->getX();
			image->modelToWorldCoordinates(coords);

			// World 2 to world 1 transformation
			coords -= origOrigin;
			coords += *newOrigin;

			wrld->worldToModelCoordinates(coords);
			atom->setX(coords);
		}
	}
	return true;
}

#ifdef BINARY

// These work with single figure or multi-object (single figure) objects.  They are all
// basically 'pass-down' functions that call the actual function in the spec'd figure.

double M3DObject::dist2FromObject(M3DObject * object, int figureID,
	DistanceType DT, bool verbose)
{
	M3DFigure * this_figure = getFigurePtr(figureID);
    M3DFigure * that_figure = object->getFigurePtr(figureID);

	double d = this_figure->dist2FromFigure(that_figure, DT, verbose);
	if (verbose) {
		if (DT == GEODESIC_DIST)             cout << "Geodesic ";
		else if (DT == EUCLIDEAN_DIST)       cout << "Euclidean ";
		else if (DT == AVE_GEODESIC_DIST)    cout << "Ave Geodesic ";
		else if (DT == AVE_EUCLIDEAN_DIST)   cout << "Ave Eudlicean ";
		cout << "Object-To-Object Distance: " << d << endl;
	}
	return d;
}

M3DObject * M3DObject::fromAveOfNeighbors(PrimNeighborhoodDefn PND, int figureID)
{
	M3DObject * retObj = assign();
	M3DFigure * this_figure = retObj->getFigurePtr(figureID);
	M3DFigure * aveFigure = this_figure->fromAveOfNeighbors(PND);
	retObj->setFigurePtr(figureID, aveFigure);
	return retObj;
}

double M3DObject::dist2FromAveOfNeighbors(PrimNeighborhoodDefn PND, int figureID,
	DistanceType DT, bool verbose)
{
	M3DObject * aveObject = fromAveOfNeighbors(PND, figureID);
	double dist2 = dist2FromObject(aveObject,figureID, DT,verbose);
	delete aveObject;
	return dist2;
}

M3DObject * M3DObject::fractionalStepToObject(M3DObject* targets[], int weights[],
	int numTargets, DistanceType DT, int figureID)
{
	M3DObject * retObj = clone();
	M3DFigure * this_figure = retObj->getFigurePtr(figureID);
	M3DFigure ** target_figs = new M3DFigure *[numTargets];

	// Create array of figures
	for (int i = 0; i < numTargets; i++)
	{
	target_figs[i] = targets[i]->getFigurePtr(figureID);
	}
 	M3DFigure * aveFigure = this_figure->fractionalStepToFigure(target_figs, weights, numTargets, DT);
	retObj->setFigurePtr(figureID,aveFigure);
	return retObj;
}

#endif

// Zhiyuan Liu
void M3DObject::renderLinkingStructure(bool isAfterSelection)
{
    int numFigs = this->getFigureCount();
    for(int i = 0; i < numFigs; ++i)
    {
        M3DQuadFigure* figure = dynamic_cast<M3DQuadFigure*>(this->getFigurePtr(i));
        std::vector<M3DSpoke> upSpokes, downSpokes;

        // 1. show interpolate 2d Srep
        figure->interpolate2D(upSpokes, 0);
        figure->interpolate2D(downSpokes, 1);

        // figure->draw2D(upSpokes, i);
        // figure->draw2D(downSpokes, i);
    }

/***************    // Computing linking structure**********************/
/***********************************************************************/
//    int numFigs = this->getFigureCount();

    // all spokes for this configuration ( 3 regions in this experiment)
    std::vector<std::vector<M3DSpoke> > configuration;
    std::vector<std::vector<M3DSpoke> > finerUpConfig;
    std::vector<std::vector<M3DSpoke> > finerDownConfig;
    std::vector<M3DQuadFigure> figures;
    std::vector<std::vector<M3DSpoke> > linkingStructure;
    for(int i = 0; i < numFigs; ++i)
    {
        M3DQuadFigure* figure = dynamic_cast<M3DQuadFigure*>(this->getFigurePtr(i));
        std::vector<M3DSpoke> upSpokes, downSpokes;

        // 1. interpolate 2d Srep in default coarse level: 3
        figure->interpolate2D(upSpokes, 0);
        figure->interpolate2D(downSpokes, 1);

        // the upSpokes and downSpokes here include two crest spokes at each vector
        int numSpokes = upSpokes.size();
        for(int j = 1; j < numSpokes - 1; ++j)
        {
            upSpokes[j].figureID = i;
            downSpokes[j].figureID = i;

            upSpokes[j].primitiveIndex = j-1;
            downSpokes[j].primitiveIndex = j-1;

            upSpokes[j].spokeType = 0;
            downSpokes[j].spokeType = 1;
        }

        // crest spokes setting
        upSpokes[0].spokeType = 2;
        upSpokes[0].primitiveIndex = 0;
        upSpokes[0].figureID = i;

        upSpokes[numSpokes-1].spokeType = 2;
        upSpokes[numSpokes-1].primitiveIndex = numSpokes-3;
        upSpokes[numSpokes-1].figureID = i;

        downSpokes[0].spokeType = 2;
        downSpokes[0].primitiveIndex = 0;
        downSpokes[0].figureID = i;

        downSpokes[numSpokes-1].spokeType = 2;
        downSpokes[numSpokes-1].figureID = i;
        downSpokes[numSpokes-1].primitiveIndex = numSpokes-3;

        // convert spokes to figure
        M3DQuadFigure tempFigure(upSpokes.size() - 2, 1);
        convertSpokes2Figure(upSpokes, downSpokes, tempFigure);
        figures.push_back(tempFigure);
        // now upSpokes contains all spokes of this figure
        upSpokes.insert(upSpokes.end(), downSpokes.begin() + 1, downSpokes.end() - 1);

        // for(int j = 0; j < upSpokes.size(); ++j)
        // {
        //     upSpokes[j].draw();
        // }
        configuration.push_back(upSpokes);

        // interpolate to finer scale 
        std::vector<M3DSpoke> finerUpSpokes, finerDownSpokes;
        int interpolateLevel = 10;
        M3DInterpolater *interp = new M3DInterpolater(figure);
        interp->interpolate2D(figure, interpolateLevel, 0, finerUpSpokes);
        interp->interpolate2D(figure, interpolateLevel, 1, finerDownSpokes);

        for(int j = 1; j < finerUpSpokes.size()-1; ++j)
        {
            finerUpSpokes[j].spokeType = 0;
            finerDownSpokes[j].spokeType = 1;
        }
        finerUpConfig.push_back(finerUpSpokes);
        finerDownConfig.push_back(finerDownSpokes);

    }

    // looking for linking correspondence for each spoke
    for(int i = 0; i < configuration.size(); ++i)
    {
        int middleMain = 1; // select middle object as main object
        int rightMain = 2;
        int leftMain = 0;
        if(i != middleMain)
            continue;
        std::vector<M3DSpoke> homeFig = configuration[i];
        for(int j = 0; j < configuration.size(); ++j)
        {
            if(j == i) continue;

            std::vector<M3DSpoke> refFig = configuration[j];
            M3DQuadFigure homeFigure, refFigure;

            homeFigure = figures[i];
            refFigure = figures[j];

            // iterate all spokes in home figure
            for(int m = 0; m < homeFig.size(); ++m)
            {
                if(m < homeFig.size()/2 && m > 0 && j== 2)
                    continue; // the left spokes in the middle main object, then don't search right object
                if(m >homeFig.size()/2-1 && m < homeFig.size() - 1 && j == 0) continue;
                M3DSpoke targetS = homeFig[m];
                Vector3D targetB = targetS.getB();
                Vector3D targetU = targetS.getU();
                bool isLinked = false;
                double delta_min = 10000.0;
                double delta_max = -1.0;
                double min_n = 0.0;
                double max_n = 0.0;
                double t_min;
                int upSpokesSize = homeFigure.getPrimitiveCount() + 2;
                int primHomeIndex = (m % upSpokesSize) - 1; // index of primitive

                // looking for the min delta from coarsly interpolated spokes
                for(int n = 0; n < refFig.size(); ++n)
                {
                    M3DSpoke ref = refFig[n];
                    Vector3D diff = ref.getB()-targetB;

                    int primRefIndex = (n % upSpokesSize) - 1; // index of primitive
                    // if unlinked
                    if(!isLinkPossible(ref, targetS)) continue;

                    isLinked = true;
                    // if the boundary overlap, then they are linked
                    if(diff.norm() < 1e-6)
                    {
                        delta_min = 0.0;
                        t_min = 0.0;
                    }
                    else
                    {
                        // prepare max and min delta for secant method
                        double t = diff*targetU / (1-ref.getU() * targetU);
                        assert(t > 0);

                        Vector3D delta_v = diff + t * (ref.getU() - targetU);

                        double delta_n = delta_v.norm();
                        double temp = delta_n * delta_n;
                        if(temp < delta_min)
                        {
                            delta_min = temp;
                            min_n = primRefIndex;
                        }
                        else if(delta_max < temp)
                        {
                            delta_max = temp;
                            max_n = primRefIndex;
                        }

                        t_min = t;
                    }

                } // end of interation on ref spokes

                // this spoke of home figure is unlinked or link to boundary
                if(!isLinked)
                {
                    continue;
                }

                /********************Now this spoke of home figure is linked to some spoke in refFig********************/
                /********************But need more work to find which exact spoke is linked to ********************/

                // already found the correspondence, no need to interpolate.
                // A:first possible outlet: this spoke of home figure is linked with existing spokes in refFigure
                // if(abs(delta_min) < 1e-6)
                // {
                //     // set linking structure for both overlap cases and meeting with existing spokes cases
                //     homeFig[m].linkTo = j;
                //     homeFig[m].linkLength = homeFig[m].getR() + t_min;
                //     homeFig[m].linkPoint = homeFig[m].getX() + homeFig[m].linkLength * targetU;
                //     homeFig[m].linkPointOnSkeleton = refFig[min_n].getX();
                //     homeFig[m].linkPointOnBdry = refFig[min_n].getB();
                // }
                // else
                {
                    // find linked spoke in finer scale
                    double delta_min = 10;
                    M3DSpoke linkedSpoke;
                    std::vector<M3DSpoke> upSpokes = finerUpConfig[j];
                    std::vector<M3DSpoke> downSpokes = finerDownConfig[j];
                    int linkSpokeType = -1;

                    for(int iS = 0; iS < upSpokes.size() - 1; ++iS)
                    {
                        Vector3D dist = upSpokes[iS].getX() - targetS.getB(); // vector from skeleton to another figure's boundary
                        double delta_temp = computeDelta(upSpokes[iS],targetS);
                        double t1 = downSpokes[iS].getU() * dist / (dist.norm()); // tester is unit vector, this is cos<tester, dist>
                        double t2 = upSpokes[iS].getU() * dist / dist.norm(); // cos<s0, dist>

                        if(delta_temp < delta_min && delta_temp >= 0)// && t1 < t2)
                        {
                            delta_min = delta_temp;
                            linkedSpoke = upSpokes[iS];
                            linkSpokeType = 0;
                            linkedSpoke.spokeId = iS;
                        }
                    }

                    for(int iS = 0; iS < downSpokes.size() - 1; ++iS)
                    {
                        Vector3D dist = downSpokes[iS].getX() - targetS.getB(); // vector from skeleton to another figure's boundary
                        double t1 = downSpokes[iS].getU() * dist / (dist.norm()); // tester is unit vector, this is cos<tester, dist>
                        double t2 = upSpokes[iS].getU() * dist / dist.norm(); // cos<s0, dist>

                        double delta_temp = computeDelta(downSpokes[iS],targetS);
                        if(delta_temp < delta_min && delta_temp >= 0)// && t1 > t2)
                        {
                            delta_min = delta_temp;
                            linkedSpoke = downSpokes[iS];
                            linkSpokeType = 1;
                            linkedSpoke.spokeId = iS;
                        }
                    }
                    double t_min = (linkedSpoke.getB() - targetB) * targetU / (1-linkedSpoke.getU() * targetU);
                    if(t_min <= 0)
                    {
                        continue;
                    }

                    double l = targetS.getR() + t_min;
                    double totalLength = l + delta_min;
                    if(delta_min / homeFig[m].linkMinDelta>1e7)
                    {
                        std::cout << "the offset greater than 1e-4, m=" << m << ", j="<<j<< std::endl;
                    }
                    if(linkSpokeType != -1 && totalLength < homeFig[m].linkLength + homeFig[m].linkMinDelta
                            && delta_min / homeFig[m].linkMinDelta<1e7)
                    {
                        homeFig[m].linkLength = l;
                        homeFig[m].linkTo = j;
                        homeFig[m].linkPoint = homeFig[m].getX() + l * targetU;
                        homeFig[m].linkPointOnSkeleton = linkedSpoke.getX();
                        homeFig[m].linkPointOnBdry = linkedSpoke.getB();
                        homeFig[m].linkSpokeType = linkSpokeType;
                        homeFig[m].isValid = true;

                        homeFig[m].linkPointOnAxis = linkedSpoke.getX() + linkedSpoke.getU() * (linkedSpoke.getR() + t_min);
                        homeFig[m].linkId = linkedSpoke.spokeId;
                        homeFig[m].linkMinDelta = delta_min;
                        
                    }

                }

            } // end of iterate in all spokes of home figure 

            // check the validity of each spoke of home figure
            Vector3D crest_home = homeFig[0].getU();
            Vector3D crest_ref = refFig[0].getU();
            int another_crest_index = refFig.size() / 2;
            Vector3D crest_ref_2 = refFig[another_crest_index].getU();

            for(int m = 1; m < homeFig.size() - 1; ++m)
            {
                if(homeFig[m].linkTo == -1) continue;
                if(homeFig[m].spokeType != 2 && homeFig[m].linkSpokeType ==  homeFig[m].spokeType && homeFig[m].spokeType != -1)
                {
                    // up spoke link to up spoke or down link to down is illegal
                    homeFig[m].isValid = false;
                }

                // if the linked spoke is crest spoke, this link is invalid also
                Vector3D linkedSpoke = homeFig[m].linkPointOnBdry - homeFig[m].linkPointOnSkeleton;
                linkedSpoke.normalize();
                if(abs(linkedSpoke * crest_ref - 1) < 1e-6 || abs(linkedSpoke * crest_ref_2 -1) < 1e-6)
                {
                    homeFig[m].isValid = false;
                }
            }

        } // end of iteration on changing secondary figure j
        
        // save new linking structure, push all objects in this configuration to linking structure
        //linkingStructure.push_back(homeFig);

        // render linking structure for home figure
        
        if(i == middleMain){
            linkingStructure.push_back(homeFig); // only select the middle object as main object

            if(selectedLinks.empty())
                drawAllSpokes(homeFig);
            else
            {
                std::vector<M3DSpoke> tempFig;
                for(int sid = 0; sid < selectedLinks.size(); ++sid)
                {
                    int linkID = selectedLinks[sid];
                    tempFig.push_back(homeFig[linkID]);
                }
                drawAllSpokes(tempFig);
            }
        }

    } // end iteration on changing primary figure i

    saveLinkingStructure(linkingStructure);
}

void M3DObject::saveLinkingStructure(std::vector<std::vector<M3DSpoke> > &linkingStructure)
{
    char filename[200];
    char *f_name = NULL;
    if(selectedLinks.empty())
    {
        f_name =  "\/diseased_linking_structure_test.m3d";
    }
    else
    {
        f_name = "\/final_links.m3d";
    }
    if(outputLinkingStructurePath != NULL)
    {
        strcpy(filename, outputLinkingStructurePath);
        strcat(filename, f_name);
    }
    else
    {
        char *default_file_path = "/playpen/workspace/demo_multiobject";
        strcpy(filename, default_file_path);
        strcat(filename, f_name);
    }

    // Sorry for the confusion of redudency.
    // It's just a terrible extension.
    if(!outputLinkingFile.empty())
    {
        strncpy(filename, outputLinkingFile.c_str(), 200);
    }
    std::cout << "Saving linking structure to file:" << filename << std::endl;
    std::ofstream f;
    f.open(filename);
    if(!f)
    {
        std::cout << "[M3DObject:saveLinkingStructure] error: can't open file:" << filename << std::endl;
        assert(true);
    }
    f << "modle {\n";
    f << "\tnumOfObj = " << linkingStructure.size() << ";\n";
    f << "\tname = default;\n";
    // f << "\tfigureTrees {\n";
    // f << "\t\tcount = 1;\n";
    // f << "\t\ttree[0]"
    for(int i = 0; i < linkingStructure.size(); ++i)
    {
        f << "\tfigure[" << i << "] {\n";
        f << "\t\tname = test;\n";

        // save links for each object
        std::vector<M3DSpoke> singleObj = linkingStructure[i];
        f << "\t\tnumOfLinks = " << singleObj.size() << ";\n";

        // save information for each single link
        for(int j = 0; j < singleObj.size(); ++j)
        {
            M3DSpoke thisLink = singleObj[j];
            f << "\t\tlink[" << j << "] {\n";

            // position
            Vector3D p = thisLink.getX();
            f << "\t\t\tx = " << p.getX() << ";\n";
            f << "\t\t\ty = " << p.getY() << ";\n";
            f << "\t\t\tz = " << p.getZ() << ";\n";

            // direction
            Vector3D u = thisLink.getU();
            f << "\t\t\tux = " << u.getX() << ";\n";
            f << "\t\t\tuy = " << u.getY() << ";\n";
            f << "\t\t\tuz = " << u.getZ() << ";\n";

            // spoke length
            f << "\t\t\tr = " << thisLink.getR() << ";\n";

            // link length from skeleton to link point
            f << "\t\t\tlinkLength = " << thisLink.linkLength << ";\n";
            int linkTo = thisLink.linkTo;
            if(!thisLink.isValid) linkTo = -1;
            f << "\t\t\tlinkTo = " << linkTo << ";\n";

            // link point z
            Vector3D z = thisLink.linkPoint;
            f << "\t\t\tlinkPoint_x = " << z.getX() << ";\n";
            f << "\t\t\tlinkPoint_y = " << z.getY() << ";\n";
            f << "\t\t\tlinkPoint_z = " << z.getZ() << ";\n";

            // spoke type: 0-up spoke 1-down spoke 2-crest spoke
            f << "\t\t\tspokeType = " << thisLink.spokeType << ";\n";

            // the vector from z' to z
            Vector3D z_prime = thisLink.linkPointOnAxis;
            Vector3D linkVector = z - z_prime;
            f << "\t\t\tlinkVector_x = " << linkVector.getX() << ";\n";
            f << "\t\t\tlinkVector_y = " << linkVector.getY() << ";\n";
            f << "\t\t\tlinkVector_z = " << linkVector.getZ() << ";\n";
            
            // linked point on the other object than the target object
            Vector3D linkedRegion = thisLink.linkPointOnBdry;
            f << "\t\t\tlinkedRegion_x = " << linkedRegion.getX() << ";\n";
            f << "\t\t\tlinkedRegion_y = " << linkedRegion.getY() << ";\n";
            f << "\t\t\tlinkedRegion_z = " << linkedRegion.getZ() << ";\n";

            f << "\t\t}\n";
        }
        f << "\t}\n";
    }

    f << "}\n";
}

void M3DObject::convertSpokes2Figure(std::vector<M3DSpoke> &upSpokes, std::vector<M3DSpoke> &downSpokes, M3DQuadFigure& output)
{
    // 1. generate a figure with input spokes
    // 1.1 first primitive is made of first 3 spokes
    M3DQuadEndPrimitive* startPrimitive = new M3DQuadEndPrimitive(upSpokes[0].getX(),   // x
                                                                upSpokes[1].getR(),   // R0
                                                                downSpokes[1].getR(), // R1
                                                                upSpokes[0].getR(),   // REnd
                                                                upSpokes[1].getU(),   // U0
                                                                downSpokes[1].getU(), // U1
                                                                upSpokes[0].getU());  // UEnd
    output.setPrimitivePtr(0, startPrimitive);

    // 1.2 last primitive is made of last 3 spokes
    int index = upSpokes.size() - 1;
    M3DQuadEndPrimitive* endPrimitive = new M3DQuadEndPrimitive(upSpokes[index].getX(),   // x
                                                                upSpokes[index-1].getR(),   // R0
                                                                downSpokes[index-1].getR(), // R1
                                                                upSpokes[index].getR(),   // REnd
                                                                upSpokes[index-1].getU(),   // U0
                                                                downSpokes[index-1].getU(), // U1
                                                                upSpokes[index].getU());  // UEnd
    output.setPrimitivePtr(output.getPrimitiveCount()-1, endPrimitive);

    // 1.3 interior primitives are standard primitives made up of 2 adjacent spokes
    for(int i = 1; i < output.getPrimitiveCount() - 1; ++i)
    {
        M3DQuadPrimitive* p = new M3DQuadPrimitive(upSpokes[i+1].getX(),     // X
                                                   upSpokes[i+1].getR(),     // R0
                                                   downSpokes[i+1].getR(),   // R1
                                                   0.0,                      // REnd
                                                   upSpokes[i+1].getU(),     // U0
                                                   downSpokes[i+1].getU(),   // U1
                                                   Vector3D(.0,.0,.0)        // UEnd

            );
        output.setPrimitivePtr(i, p);
    }

}
// in parameter space, use secant method to find 0-crossing of delta
// the second point is (min_n, delta_min), the first point is (max_n, delta_max)
// given a position x on 2D s-rep, need to output delta and that spoke:
// if x is integer, then just use original spokes,
// otherwise, need interpolate.
double M3DObject::computeDelta2D(M3DFigure& fig, M3DSpoke& targetS, double x, M3DSpoke& output)
{
    // 2. find x in the parameter space of new figure above

    // 2.1 if x is integer, then use original spoke as y'
    if(abs(x-int(x)) < 1e-6)
    {
        M3DQuadPrimitive* p = dynamic_cast<M3DQuadPrimitive*>(fig.getPrimitivePtr(x));
        output.setX(p->getX());
        M3DSpoke s0(p->getX(), p->getU0(), p->getR0());
        M3DSpoke s1(p->getX(), p->getU1(), p->getR1());
        return findBetterLinkedSpoke(targetS, s0, s1, output);
    }

    // 2.2 if x is not integer, interpolate fig first, then compute delta
    M3DInterpolater *interp = new M3DInterpolater(&fig);
    M3DSpoke upSpoke, downSpoke;
    interp->getSpokeFromInterpolation(&fig, 0, x, upSpoke);
    interp->getSpokeFromInterpolation(&fig, 1, x, downSpoke);
    return findBetterLinkedSpoke(targetS, upSpoke, downSpoke, output);

}
double M3DObject::findBetterLinkedSpoke(M3DSpoke& targetS, M3DSpoke& s0, M3DSpoke& s1, M3DSpoke& output)
{
    double delta0 = computeDelta(s0, targetS);
    double delta1 = computeDelta(s1, targetS);

    if(delta0 < 0 && delta1 >= 0)
    {
        output.copyFrom(s1);
        return delta1;
    }
    else if(delta0 < 0 && delta1 < 0)
    {
        return -1.0;
    }
    else if(delta0 >= 0 && delta1 >= 0)
    {
        if(delta0 < delta1)
        {
            output.copyFrom(s0);
            return delta0;
        }
        else
        {
            output.copyFrom(s1);
            return delta1;
        }

    }
    else if(delta0 >= 0 && delta1 < 0)
    {
        output.copyFrom(s0);
        return delta0;
    }

}

// fixed y_prime, varing y in neighboring objects
double M3DObject::computeDelta(M3DSpoke& y_prime, M3DSpoke& y)
{
    Vector3D b_prime = y_prime.getB();
    Vector3D b = y.getB();
    Vector3D u_prime = y_prime.getU();
    Vector3D u = y.getU();
    double t = (b_prime - b) * u / (1-u_prime * u);
    if(t < 0 || !isLinkPossible(y_prime, y))
    {
        return -1.0;
    }
    double delta = ((b_prime - b) + t * (u_prime - u)).norm();
    return delta * delta;
}
bool M3DObject::isLinkPossible(M3DSpoke& s0, M3DSpoke& s1)
{
    Vector3D y_prime = s0.getB();
    Vector3D y = s1.getB();
    Vector3D diff = y_prime - y; // vector from bdry to bdry
    Vector3D dist = s0.getX() - s1.getB(); // vector from skeleton to another figure's boundary
    double dist_norm_sq = dist.norm() * dist.norm();

    if(diff * s1.getU() > 0               //(y'-y) * U > 0
       && (diff * dist > 0 && diff * dist < dist_norm_sq)// vector y' - y shouldnot pass through ref skeleton
        )
    {
        return true;
    }
    return false;

}

void M3DObject::drawAllSpokes(std::vector<M3DSpoke> &homeFig)
{
    bool showThese = true;

    // draw external linking axis
    for(int m = 0; m < homeFig.size(); ++m)
    {
        //showThese = ((m >= 19 && m <= 48));
        if(homeFig[m].linkTo == -1) continue;

        if(!showThese) continue;
        double rgb[3] = {0.0, 1.0, 0.0};
        if(!homeFig[m].isValid)
        {
            // invalid display as red
            rgb[0] = 1.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;
        }
        else
        {

        }
        Vector3D axisPoint = homeFig[m].linkPoint;
        axisPoint.draw(0.005, rgb[0], rgb[1], rgb[2]);
    }

    // draw connection from linking axis to home figure
    glColor3f(0.5f, 0.5f, 1.0f);
    glBegin(GL_LINES);

    for(int m = 0; m < homeFig.size(); ++m)
    {
        //showThese = ((m >= 19 && m <= 48));
        if(homeFig[m].linkTo == -1 || !showThese) continue;
        Vector3D axisPoint = homeFig[m].linkPoint;
        glVertex3f(axisPoint.getX(), axisPoint.getY(), axisPoint.getZ());
        Vector3D medialPoint = homeFig[m].getX();
        glVertex3f(medialPoint.getX(), medialPoint.getY(), medialPoint.getZ());

    }
    glEnd();

    // draw spokes on home figure
//    glLineWidth(2.0);
//    glBegin(GL_LINES);

//    for(int m = 0; m < homeFig.size(); ++m)
//    {
//        if(homeFig[m].linkTo == -1) continue;

//        if(m == 34)
//        {
//           glColor3f(0.0f, 0.2f, 0.5f);

//        }
//        Vector3D impliedBdry = homeFig[m].getB();
//        glVertex3f(impliedBdry.getX(), impliedBdry.getY(), impliedBdry.getZ());
//        Vector3D medialPoint = homeFig[m].getX();
//        glVertex3f(medialPoint.getX(), medialPoint.getY(), medialPoint.getZ());

//    }
//    glEnd();

    // draw connection from linking axis to reference figure
    glColor3f(0.5f, 0.5f, 1.0f);
    glBegin(GL_LINES);

    for(int m = 0; m < homeFig.size(); ++m)
    {
        //showThese = ((m >= 19 && m <= 48));
        if(homeFig[m].linkTo == -1 || !showThese) continue;

        Vector3D axisPoint = homeFig[m].linkPoint;
        glVertex3f(axisPoint.getX(), axisPoint.getY(), axisPoint.getZ());
        Vector3D linkedPoint = homeFig[m].linkPointOnSkeleton;
        glVertex3f(linkedPoint.getX(), linkedPoint.getY(), linkedPoint.getZ());

    }
    glEnd();

    // // draw the linked spoke in ref figure
//     glColor3f(.0f, 1.0f, 0.0f);
//     glBegin(GL_LINES);


//     for(int m = 0; m < homeFig.size(); ++m)
//     {
//         if(homeFig[m].linkTo == -1 || (m != 39 && m != 40)) continue;

//         if(m == 4) glColor3f(1.0f, 1.0f, 1.0f);
//         if(m==5) glColor3f(.0f, 1.0f, 1.0f);
//         Vector3D bdryPoint = homeFig[m].linkPointOnBdry;
//         glVertex3f(bdryPoint.getX(), bdryPoint.getY(), bdryPoint.getZ());
//         Vector3D linkedPoint = homeFig[m].linkPointOnSkeleton;
//         glVertex3f(linkedPoint.getX(), linkedPoint.getY(), linkedPoint.getZ());
//     }
//     glEnd();


}

void M3DObject::readSelectionResult(const std::string& filename)
{
    std::vector<int> idList;
    std::ifstream infile(filename.c_str());
    int id;
    while(infile >> id)
    {
        idList.push_back(id);
    }

    if(this->selectedLinks.empty())
    {
        this->selectedLinks = idList;
    }
    else
    {
        this->selectedLinks.clear();
    }
}
