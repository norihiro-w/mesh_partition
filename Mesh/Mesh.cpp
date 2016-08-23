#include "Mesh.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "Node.h"
#include "Edge.h"
#include "Elem.h"

//#define BUILD_MESH_EDGE

using namespace std;

namespace
{

template<typename T> std::string number2str(T d)
{
	std::stringstream out;
	out << d;
	return out.str();
}

void read_npart_file(const int num_parts, const string fname, vector<long> &vec_node_dom_idx, vector<bool> &vec_node_dom_marked)
{
	const string s_nparts(number2str(num_parts));
	const string f_iparts = fname + ".mesh.npart." + s_nparts;
	ifstream npart_in(f_iparts.c_str());
	if (!npart_in.is_open())
	{
		cerr << ("Error: cannot open .npart file . It may not exist !");
		exit(1);
	}
	for (long i = 0; i < static_cast<long>(vec_node_dom_idx.size()); i++)
	{
		int dom;
		npart_in >> dom >> ws;
		vec_node_dom_idx[i] = dom;
		vec_node_dom_marked[i] = false;
	}
	npart_in.close();
}

} // namespace


//------------------------------------------------------
//   Topology definition of geometrical element.
//    WW. 10.01.2005
//------------------------------------------------------
namespace Mesh_Group
{

Mesh::Mesh(bool quad)
{
	useQuadratic = quad;
	coordinate_system = 1;
	axisymmetry = false;
	max_ele_dim = 0;
	NodesNumber_Linear = NodesNumber_Quadratic = 0;
}

Mesh::~Mesh()
{
	long i;
	// Nodes
	for (i = 0; i < (long) node_vector.size(); i++)
		delete node_vector[i];
	node_vector.clear();
	// Edges
#ifdef BUILD_MESH_EDGE
	for(i=0; i<(long)edge_vector.size(); i++)
	delete edge_vector[i];
	edge_vector.clear();
#endif

#ifdef BUILD_MESH_FACE
	// Surface faces
	for(i=0; i<(long)face_vector.size(); i++)
	delete face_vector[i];
	face_vector.clear();
#endif

	// Element
	for (i = 0; i < (long) elem_vector.size(); i++)
		delete elem_vector[i];
	elem_vector.clear();
}

// Construct grid
//
/**************************************************************************
 ConnectedNodes
 **************************************************************************/
void Mesh::ConnectedNodes(bool quadratic)
{
	int i, j, l, k, n;
	Node* m_nod = NULL;
	Elem* m_ele = NULL;
	bool exist = false;
	//----------------------------------------------------------------------
	for (i = 0; i < (long) node_vector.size(); i++)
	{
		m_nod = node_vector[i];
		for (j = 0; j < (int) m_nod->ElementsRelated.size(); j++)
		{
			m_ele = elem_vector[m_nod->ElementsRelated[j]];
			for (l = 0; l < m_ele->getNodesNumber(quadratic); l++)
			{
				exist = false;
				const long nidx = m_ele->nodes[l]->index;
				for (k = 0; k < (int) m_nod->NodesRelated.size(); k++)
				{
					if (m_nod->NodesRelated[k] == nidx)
					{
						exist = true;
						break;
					}
				}
				if (!exist)
					m_nod->NodesRelated.push_back(nidx);
			}
		}
	}
	for (i = 0; i < (long) node_vector.size(); i++)
	{
		m_nod = node_vector[i];
		j = (int) m_nod->NodesRelated.size();
		for (k = 0; k < j; k++)
		{
			for (l = k; l < j; l++)
			{
				if (m_nod->NodesRelated[l] < m_nod->NodesRelated[k])
				{
					n = m_nod->NodesRelated[k];
					m_nod->NodesRelated[k] = m_nod->NodesRelated[l];
					m_nod->NodesRelated[l] = n;
				}
			}
		}

	}
}

/**************************************************************************
 ConnectedElements2Node
 **************************************************************************/
void Mesh::ConnectedElements2Node(bool quadratic)
{
	for (auto node :  node_vector)
		node->ElementsRelated.clear();
	// set neighbors of node
	for (long e = 0; e < (long) elem_vector.size(); e++)
	{
		auto thisElem0 = elem_vector[e];
		if (!thisElem0->getStatus())
			continue;      // Not marked for use
		for (int i = 0; i < thisElem0->getNodesNumber(quadratic); i++)
		{
			bool found = false;
			auto node = node_vector[thisElem0->getNodeIndex(i)];
			for (int j = 0; j < (int) node->ElementsRelated.size(); j++)
			{
				if (e == node->ElementsRelated[j])
				{
					found = true;
					break;
				}
			}
			if (!found)
				node->ElementsRelated.push_back(e);
		}
	}
}

void Mesh::ConstructGrid()
{
	int faceIndex_loc0[10];
	int faceIndex_loc[10];
	std::vector<Node*> e_nodes0(20);
	long node_index_glb[20];
	long node_index_glb0[20];

#ifdef BUILD_MESH_EDGE
	int edgeIndex_loc0[3];
	int edgeIndex_loc[3];
	std::vector<int> Edge_Orientation(15, 1);
	std::vector<Edge*> Edges(15);
	std::vector<Edge*> Edges0(15);
#endif
	std::vector<Elem*> Neighbors(15);
	std::vector<Elem*> Neighbors0(15);

	std::vector<Node*> e_edgeNodes0(3);
	std::vector<Node*> e_edgeNodes(3);

	clock_t start, finish;
	start = clock();

	//Elem->nodes not initialized

	if (NodesNumber_Linear == 0)
		NodesNumber_Linear = (long) node_vector.size();

	//----------------------------------------------------------------------
	// set neighbors of node
	ConnectedElements2Node(useQuadratic);
	//----------------------------------------------------------------------

	//----------------------------------------------------------------------
	// Compute neighbors and edges
	auto const n_elements = (long) elem_vector.size();
	for (long e0_id = 0; e0_id < n_elements; e0_id++)
	{
		auto thisElem0 = elem_vector[e0_id];
		thisElem0->setOrder(useQuadratic);
		// get nodes
		const int nnodes0 = thisElem0->getNodesNumber(useQuadratic);
		thisElem0->getNodeIndeces(node_index_glb0);
		for (int i = 0; i < nnodes0; i++)
			e_nodes0[i] = node_vector[node_index_glb0[i]];
		// get neighbors
		thisElem0->getNeighbors(Neighbors0);
		const int nElem0Faces = thisElem0->getFacesNumber();
		// set neighbors
		for (int i = 0; i < nElem0Faces; i++)
		{
			if (Neighbors0[i])
				continue;

			// look for an element sharing the same face
			bool foundNeighbor = false;
			const int nElem0FaceNodes = thisElem0->getElementFaceNodes(i, faceIndex_loc0);
			for (int k = 0; k < nElem0FaceNodes; k++)
			{
				Mesh_Group::Node* face_node = e_nodes0[faceIndex_loc0[k]];
				const long n_elems_connected_to_face_node = (long) face_node->ElementsRelated.size();
				for (long ei = 0; ei < n_elems_connected_to_face_node; ei++)
				{
					const long conn_ele_id = face_node->ElementsRelated[ei];
					if (conn_ele_id == e0_id)
						continue; //skip same element

					Mesh_Group::Elem* connectedElem = elem_vector[conn_ele_id];
					connectedElem->getNodeIndeces(node_index_glb);
					connectedElem->getNeighbors(Neighbors);
					const int nConnElemFaces = connectedElem->getFacesNumber();

					for (int ii = 0; ii < nConnElemFaces; ii++) // Faces of the connected element
					{
						const int nConnElemFaceNodes = connectedElem->getElementFaceNodes(ii, faceIndex_loc);
						// check if this face is shared
						if (nElem0FaceNodes != nConnElemFaceNodes)
							continue;
						int counter = 0;
						for (int j = 0; j < nElem0FaceNodes; j++)
						{
							for (int jj = 0; jj < nConnElemFaceNodes; jj++)
							{
								if (node_index_glb0[faceIndex_loc0[j]] == node_index_glb[faceIndex_loc[jj]])
								{
									counter++;
									break;
								}
							}
						}
						if (counter != nConnElemFaceNodes)
							continue;
						// found neighbor for this face
						Neighbors0[i] = connectedElem;
						Neighbors[ii] = thisElem0;
						connectedElem->setNeighbor(ii, thisElem0);
						foundNeighbor = true;
					}
					if (foundNeighbor)
						break;
				}
				if (foundNeighbor)
					break;
			}
		}
		thisElem0->setNeighbors(Neighbors0);

#ifdef BUILD_MESH_EDGE
		// --------------------------------
		// Edges
		const int nedges0 = thisElem0->getEdgesNumber();
		thisElem0->getEdges(Edges0);
		for(int i=0; i<nedges0; i++)
		{
			thisElem0->getLocalIndices_EdgeNodes(i, edgeIndex_loc0);
			// Check neighbors
			bool done = false;
			for(int k=0; k<2; k++)
			{
				auto edge_node = e_nodes0[edgeIndex_loc0[k]];
				const long nConnElements = (long)edge_node->ElementsRelated.size();
				for(long ei=0; ei<nConnElements; ei++)
				{
					auto const connected_element_id = edge_node->ElementsRelated[ei];
					if(connected_element_id == e0_id) continue;
					auto connected_element = elem_vector[connected_element_id];
					connected_element->getNodeIndeces(node_index_glb);
					connected_element->getEdges(Edges);
					// Edges of neighbors
					auto const nConnEleEdges = connected_element->getEdgesNumber();
					for(int ii=0; ii<nConnEleEdges; ii++)
					{
						connected_element->getLocalIndices_EdgeNodes(ii, edgeIndex_loc);
						if(( node_index_glb0[edgeIndex_loc0[0]]==node_index_glb[edgeIndex_loc[0]]
										&&node_index_glb0[edgeIndex_loc0[1]]==node_index_glb[edgeIndex_loc[1]])
								||( node_index_glb0[edgeIndex_loc0[0]]==node_index_glb[edgeIndex_loc[1]]
										&&node_index_glb0[edgeIndex_loc0[1]]==node_index_glb[edgeIndex_loc[0]]) )
						{
							if(Edges[ii])
							{
								Edges0[i] = Edges[ii];
								Edges[ii]->getNodes(e_edgeNodes);
								if( node_index_glb0[edgeIndex_loc0[0]]==e_edgeNodes[1]->getIndex()
										&& node_index_glb0[edgeIndex_loc0[1]]==e_edgeNodes[0]->getIndex())
								Edge_Orientation[i] = -1;
								done = true;
								break;
							}
						}
					} //  for(ii=0; ii<nedges; ii++)
					if(done) break;
				} // for(ei=0; ei<e_size_l; ei++)
				if(done) break;
			} //for(k=0;k<2;k++)
			if(!done)// new edges and new node
			{
				Edges0[i] = new Edge((long)edge_vector.size());
				Edges0[i]->setOrder(false);
				e_edgeNodes0[0] = e_nodes0[edgeIndex_loc0[0]];
				e_edgeNodes0[1] = e_nodes0[edgeIndex_loc0[1]];
				e_edgeNodes0[2] = NULL;
				Edges0[i]->setNodes(e_edgeNodes0);
				edge_vector.push_back(Edges0[i]);
			} // new edges
		} //  for(i=0; i<nedges0; i++)
		  //
		  // set edges nodes
		thisElem0->setEdges_Orientation(Edge_Orientation);
		thisElem0->setEdges(Edges0);

#endif  //BUILD_MESH_EDGE
//	  // set nodes
//      thisElem0->setOrder(false);
		// Resize is true
		thisElem0->setNodes(e_nodes0, true);
	}      // Over elements

	// set faces on surfaces and others
	msh_no_line = 0;  // Should be members of mesh
	msh_no_quad = 0;
	msh_no_hexs = 0;
	msh_no_tris = 0;
	msh_no_tets = 0;
	msh_no_pris = 0;
	msh_no_pyra = 0;
	for (auto thisElem0 : elem_vector)
	{
		switch (thisElem0->getElementType())
		{
		case line:
			msh_no_line++;
			break;
		case quadri:
			msh_no_quad++;
			break;
		case hex:
			msh_no_hexs++;
			break;
		case tri:
			msh_no_tris++;
			break;
		case tet:
			msh_no_tets++;
			break;
		case prism:
			msh_no_pris++;
			break;
		case pyramid:
			msh_no_pyra++;
			break;
		}
		// Compute volume meanwhile
		//thisElem0->ComputeVolume();

		if (thisElem0->getElementType() == line)
			continue; // line element
		thisElem0->getNodeIndeces(node_index_glb0);
		thisElem0->getNeighbors(Neighbors0);
		auto const nElem0Faces = thisElem0->getFacesNumber();

#ifdef BUILD_MESH_FACE
		// Check face on surface
		for(i=0; i<m0; i++)// Faces
		{
			if(Neighbors0[i])
			continue;
			Elem* newFace = new Elem((long)face_vector.size(), thisElem0, i);
//          thisElem0->boundary_type='B';
			thisElem0->no_faces_on_surface++;
			face_vector.push_back(newFace);
			Neighbors0[i] = newFace;
		}
#endif
		thisElem0->setNeighbors(Neighbors0);

	}
	NodesNumber_Quadratic = (long) node_vector.size();
	if ((msh_no_hexs + msh_no_tets + msh_no_pris + msh_no_pyra) > 0)
		max_ele_dim = 3;
	else if ((msh_no_quad + msh_no_tris) > 0)
		max_ele_dim = 2;
	else
		max_ele_dim = 1;
	//----------------------------------------------------------------------
	// Node information
	// 1. Default node index <---> eqs index relationship
	// 2. Coordiate system flag
	double x_sum = 0.0;
	double y_sum = 0.0;
	double z_sum = 0.0;
	for (long e0_id = 0; e0_id < (long) node_vector.size(); e0_id++)
	{
		x_sum += fabs(node_vector[e0_id]->X());
		y_sum += fabs(node_vector[e0_id]->Y());
		z_sum += fabs(node_vector[e0_id]->Z());
	}
	if (x_sum > 0.0 && y_sum < DBL_MIN && z_sum < DBL_MIN)
		coordinate_system = 10;
	else if (y_sum > 0.0 && x_sum < DBL_MIN && z_sum < DBL_MIN)
		coordinate_system = 11;
	else if (z_sum > 0.0 && x_sum < DBL_MIN && y_sum < DBL_MIN)
		coordinate_system = 12;
	else if (x_sum > 0.0 && y_sum > 0.0 && z_sum < DBL_MIN)
		coordinate_system = 21;
	else if (x_sum > 0.0 && z_sum > 0.0 && y_sum < DBL_MIN)
		coordinate_system = 22;
	else if (x_sum > 0.0 && y_sum > 0.0 && z_sum > 0.0)
		coordinate_system = 32;
	/*  // 23.05.2008. WW. Futher test is needed
	 // 1D in 2D
	 if(msh_no_line>0)   //
	 {
	 if(x_sum>0.0&&y_sum>0.0&&z_sum<DBL_MIN)
	 coordinate_system = 22;
	 if(x_sum>0.0&&z_sum>0.0&&y_sum<DBL_MIN)
	 coordinate_system = 22;
	 }
	 */
	//----------------------------------------------------------------------
	// For sparse matrix
	ConnectedNodes(false);
	//
	e_nodes0.resize(0);

#ifdef BUILD_MESH_EDGE
	Edge_Orientation.resize(0);
	Edges.resize(0);
	Edges0.resize(0);
	e_edgeNodes0.resize(0);
	e_edgeNodes.resize(0);
#endif

	Neighbors.resize(0);
	Neighbors0.resize(0);

	finish = clock();
	cout << "\nCPU time elapsed in constructing topology of grids: " << (double) (finish - start) / CLOCKS_PER_SEC << "s" << endl << endl;

}

/**************************************************************************
 Programing:
 07/2007 WW Implementation
 **************************************************************************/
void Mesh::GenerateHighOrderNodes()
{
	int edgeIndex_loc0[3];
	int edgeIndex_loc1[3];


	// Set neighbors of node. All elements, even in deactivated subdomains, are taken into account here.
	for (long e = 0; e < (long) node_vector.size(); e++)
		node_vector[e]->getConnectedElementIDs().clear();
	bool done = false;
	const long ele_vector_size(elem_vector.size());
	for (long e = 0; e < ele_vector_size; e++)
	{
		Elem* thisElem0 = elem_vector[e];
		for (int i = 0; i < thisElem0->getNodesNumber(false); i++)
		{
			done = false;
			long ni = thisElem0->getNodeIndex(i);
			size_t n_connected_elements(node_vector[ni]->getConnectedElementIDs().size());
			for (size_t j = 0; j < n_connected_elements; j++)
				if (e == node_vector[ni]->getConnectedElementIDs()[j])
				{
					done = true;
					break;
				}
			if (!done)
				node_vector[ni]->getConnectedElementIDs().push_back(e);
		}
	}

	clock_t start, finish;
	start = clock();

	//
	std::vector<Node*> e_nodes0(20);
	std::vector<Node*> e_nodes(20);
	//----------------------------------------------------------------------
	NodesNumber_Linear = (long) node_vector.size();
	// Loop over elements
	const auto e_size = (long) elem_vector.size();
	bool hasLines = false;
	for (long e = 0; e < e_size; e++)
	{
		auto thisElem0 = elem_vector[e];
		if (thisElem0->getElementType() == Mesh_Group::line) {
			hasLines = true;
			continue;
		}
		auto nnodes0 = thisElem0->nnodes; // Number of nodes for linear element
		//thisElem0->GetNodeIndeces(node_index_glb0);
		for (int i = 0; i < nnodes0; i++) // Nodes
			e_nodes0[i] = thisElem0->getNode(i);
		// --------------------------------
		// Edges
		auto nedges0 = thisElem0->getEdgesNumber();
		// Check if there is any neighbor that has new middle points
		for (int i = 0; i < nedges0; i++)
		{
#ifdef BUILD_MESH_EDGE
			auto thisEdge0 = thisElem0->getEdge(i);
#endif
			thisElem0->getLocalIndices_EdgeNodes(i, edgeIndex_loc0);
			const long ena0 = thisElem0->getNodeIndex(edgeIndex_loc0[0]);
			const long ena1 = thisElem0->getNodeIndex(edgeIndex_loc0[1]);
			// Check neighbors
			bool done = false;
			for (int k = 0; k < 2; k++)
			{
				auto const nEdgeConnectedElements = (long) e_nodes0[edgeIndex_loc0[k]]->ElementsRelated.size();
				for (int ei = 0; ei < nEdgeConnectedElements; ei++)
				{
					auto ee = e_nodes0[edgeIndex_loc0[k]]->ElementsRelated[ei];
					if (ee == e)
						continue;
					auto thisElem = elem_vector[ee];
					auto nedges = thisElem->getEdgesNumber();

					// If this element already proccessed
					if (thisElem->nodes.size() == thisElem->getNodesNumberHQ())
					{
						// Edges of neighbors
						for (int ii = 0; ii < nedges; ii++)
						{
							thisElem->getLocalIndices_EdgeNodes(ii, edgeIndex_loc1);

							const long enb0 = thisElem->getNodeIndex(edgeIndex_loc1[0]);
							const long enb1 = thisElem->getNodeIndex(edgeIndex_loc1[1]);

							if (((ena0 == enb0) && (ena1 == enb1)) || ((ena0 == enb1) && (ena1 == enb0)))
							{
								auto aNode = thisElem->getNode(edgeIndex_loc1[2]);
								e_nodes0[edgeIndex_loc0[2]] = aNode;
								done = true;
								break;
							}

						} //  for(ii=0; ii<nedges; ii++)

					}

					if (done)
						break;
				} // for(ei=0; ei<e_size_l; ei++)
				if (done)
					break;
			} //for(k=0;k<2;k++)
			if (!done)
			{
				auto aNode = new Node((long) node_vector.size());
				const Node *na = thisElem0->getNode(edgeIndex_loc0[0]);
				const Node *nb = thisElem0->getNode(edgeIndex_loc0[1]);
				aNode->setX(0.5 * (na->X() + nb->X()));
				aNode->setY(0.5 * (na->Y() + nb->Y()));
				aNode->setZ(0.5 * (na->Z() + nb->Z()));
				e_nodes0[edgeIndex_loc0[2]] = aNode;

#ifdef BUILD_MESH_EDGE
				thisEdge0->setNode(2, aNode);
#endif
				node_vector.push_back(aNode);
			}
		} //  for(i=0; i<nedges0; i++)

		// No neighors or no neighbor has new middle point
		//
		if (thisElem0->getElementType() == quadri) // Quadrilateral
		{
			double x0 = 0.0, y0 = 0.0, z0 = 0.0;
			auto aNode = new Node((long) node_vector.size());
			e_nodes0[8] = aNode;
			nnodes0 = thisElem0->nnodes;
			for (int i = 0; i < nnodes0; i++) // Nodes
			{
				x0 += e_nodes0[i]->X();
				y0 += e_nodes0[i]->Y();
				z0 += e_nodes0[i]->Z();
			}
			x0 /= (double) nnodes0;
			y0 /= (double) nnodes0;
			z0 /= (double) nnodes0;
			aNode->setX(x0);
			aNode->setY(y0);
			aNode->setZ(z0);
			node_vector.push_back(aNode);
		}
		// Set edges and nodes
		thisElem0->setOrder(true);
		// Resize is true
		thisElem0->setNodes(e_nodes0, true);
	} // Over elements
	  //

	// Setup 1d line elements at the end
	if (hasLines)
	{
		for (long e = 0; e < e_size; e++)
		{
			if (elem_vector[e]->getElementType() != Mesh_Group::line)
				continue;
			auto thisEdge0 = elem_vector[e];

			const auto nnodes0 = thisEdge0->nnodes;
			for (int i = 0; i < nnodes0; i++)
				e_nodes0[i] = thisEdge0->getNode(i);

			std::vector<int> elementIDs_connected_to_edge_nodes;
			for (int i=0; i<nnodes0; i++)
				for (auto eid : node_vector[thisEdge0->getNodesNumber(i)]->ElementsRelated)
					elementIDs_connected_to_edge_nodes.push_back(eid);
			std::sort(elementIDs_connected_to_edge_nodes.begin(), elementIDs_connected_to_edge_nodes.end());
			elementIDs_connected_to_edge_nodes.erase(std::unique(elementIDs_connected_to_edge_nodes.begin(), elementIDs_connected_to_edge_nodes.end()), elementIDs_connected_to_edge_nodes.end());

			bool foundCommonEdge = false;
			for (auto ele_id : elementIDs_connected_to_edge_nodes)
			{
				auto connElem = elem_vector[ele_id];
				if (connElem->getElementType() == Mesh_Group::line)
					continue;

				// check if it has a common edge
				for (int i_edge=0; i_edge<connElem->getEdgesNumber(); i_edge++)
				{
					auto thisEdge = connElem->getEdge(i_edge);
					for (int i = 0; i < 2; i++)
						e_nodes[i] = thisEdge->getNode(i);
					bool foundEdge = true;
					for (int i=0; i<2; i++)
					{
						bool foundNode = false;
						for (int j=0; j<2; j++)
						{
							if (e_nodes[j]==e_nodes0[i])
							{
								foundNode = true;
								break;
							}
						}
						if (!foundNode) {
							foundEdge = false;
							break;
						}
					}
					if (!foundEdge)
						continue;
					//the edge is found now
					auto aNode = thisEdge->getNode(2);
					if (aNode) // The middle point exist
					{
						e_nodes0[nnodes0] = aNode;
						foundCommonEdge = true;
						break;
					}
				}
				if (foundCommonEdge)
					break;
			}

			if (!foundCommonEdge)
			{
				auto aNode = new Node((long) node_vector.size());
				double x0 = 0.0, y0 = 0.0, z0 = 0.0;
				for (int i = 0; i < nnodes0; i++) // Nodes
				{
					x0 += e_nodes0[i]->X();
					y0 += e_nodes0[i]->Y();
					z0 += e_nodes0[i]->Z();
				}
				x0 /= (double) nnodes0;
				y0 /= (double) nnodes0;
				z0 /= (double) nnodes0;
				aNode->setX(x0);
				aNode->setY(y0);
				aNode->setZ(z0);
				e_nodes0[nnodes0] = aNode;
				node_vector.push_back(aNode);
			}
			thisEdge0->setOrder(true);
			thisEdge0->setNodes(e_nodes0, true);
		}
	}

	NodesNumber_Quadratic = (long) node_vector.size();
	for (long e = 0; e < e_size; e++)
	{
		auto thisElem0 = elem_vector[e];
		for (int i = thisElem0->nnodes; i < thisElem0->nnodesHQ; i++)
		{
			bool done = false;
			auto aNode = thisElem0->getNode(i);
			for (int k = 0; k < (int) aNode->ElementsRelated.size(); k++)
			{
				if (e == aNode->ElementsRelated[k])
				{
					done = true;
					break;
				}
			}
			if (!done)
				aNode->ElementsRelated.push_back(e);
		}
	}

	// For sparse matrix
	ConnectedNodes(true);
	//ConnectedElements2Node(true);
	//
	e_nodes0.resize(0);

	finish = clock();
	cout << "\n\tCPU time elapsed in generating high oder elements: " << (double) (finish - start) / CLOCKS_PER_SEC << "s" << endl;

}

void Mesh::ConstructSubDomain_by_Elements(const string fname, const int num_parts, const bool osdom)
{
	string str;
	string stro;
	int dom;
	int max_dom;
	int k, kk;
	long i, j;
	//  int ntags = 3;

	fstream gmsh_out;

	string deli = " ";
	//

	string s_nparts;
	stringstream ss;
	ss << num_parts;
	ss >> s_nparts;
	ss.clear();

	str = fname + ".mesh.epart." + s_nparts;
	stro = fname + "." + s_nparts + "ddc";

	//namef = ".mesh.epart."; //+str_buf;
	ifstream part_in;
	fstream part_out;
	part_out.open(stro.c_str(), ios::out | ios::trunc);
	// Output for gmsh

	if (osdom)
	{
		stro = fname + "_gmsh.msh";
		gmsh_out.open(stro.c_str(), ios::out);
		//gmsh_out<<"$NOD"<<endl;
		gmsh_out << "$MeshFormat\n2 0 8\n$EndMeshFormat\n$Nodes" << endl;
		gmsh_out << node_vector.size() << endl;
		Node *node;
		for (i = 0; i < (long) node_vector.size(); i++)
		{
			gmsh_out << i + 1 << " ";
			node = node_vector[i];
			gmsh_out << node->X() << " ";
			gmsh_out << node->Y() << " ";
			gmsh_out << node->Z() << endl;
		}
		//gmsh_out<<"$ENDNOD"<<endl;
		//gmsh_out<<"$ELM"<<endl;
		gmsh_out << "$EndNodes\n$Elements" << endl;
		gmsh_out << (long) elem_vector.size() << endl;
	}

	//
	part_in.open(str.c_str());
	if (!part_in.is_open())
	{
		cerr << ("Error: cannot open .epart file . It may not exist !");
		abort();
	}

	max_dom = 0;
	Elem *ele = NULL;

	for (i = 0; i < (long) elem_vector.size(); i++)
	{
		part_in >> dom >> ws;
		ele = elem_vector[i];
		ele->setDomainIndex(dom);
//      elem_vector[i]->AllocateLocalIndexVector();
		if (dom > max_dom)
			max_dom = dom;

		if (osdom)
		{
			ele->WriteGmsh(gmsh_out, dom + 1);
		}

	}
	max_dom++;
	part_in.close();
	remove(str.c_str());

	if (osdom)
	{
		gmsh_out << "$EndElements" << endl;
		gmsh_out.close();
	}
	//

	//Output ddc file
	// long *nod_dom = new long[max_dom];
	long *ele_dom = new long[max_dom];
	for (k = 0; k < max_dom; k++)
	{
		ele_dom[k] = 0;
		//nod_dom[k]=0;
		for (j = 0; j < (long) elem_vector.size(); j++)
		{
			if (elem_vector[j]->getDomainIndex() == k)
				ele_dom[k] += 1;
		}
	}

	bool done = false;
	long n_index = 0;
	vector<int> nodes_dom;
	vector<Elem*> eles_dom;

	//
	for (k = 0; k < max_dom; k++)
	{
		part_out << "#DOMAIN " << k << endl;
		part_out << "$ELEMENTS " << ele_dom[k] << endl;
		nodes_dom.clear();
		eles_dom.clear();
		for (j = 0; j < (long) elem_vector.size(); j++)
		{
			ele = elem_vector[j];
			for (kk = 0; kk < ele->getNodesNumber(useQuadratic); kk++)
			{
				ele->setLocalNodeIndex(kk, -1);
				ele->AllocateLocalIndexVector();
				ele->setDomNodeIndex(kk, -1);
			}
		}
		for (j = 0; j < (long) elem_vector.size(); j++)
		{
			ele = elem_vector[j];
			//ele->AllocateLocalIndexVector();
			if (ele->getDomainIndex() == k)
			{
				for (kk = 0; kk < ele->getNodesNumber(useQuadratic); kk++)
				{
					done = false;
					n_index = ele->getLocalNodeIndex(kk);
					if (n_index > -1)
					{
						ele->setDomNodeIndex(kk, n_index);
						done = true;
					}
					if (!done)
					{
						ele->setDomNodeIndex(kk, (long) nodes_dom.size()); //For test output
						ele->setLocalNodeIndex(kk, (long) nodes_dom.size());
						nodes_dom.push_back(ele->getNodeIndex(kk));
					}
				}
				part_out << ele->getIndex() << endl;
				eles_dom.push_back(ele); //TEST OUT
			}
		}
		part_out << "$NODES_INNER " << (long) nodes_dom.size() << endl;
		for (j = 0; j < (long) nodes_dom.size(); j++)
			part_out << nodes_dom[j] << endl;

		if (osdom)
		{
			string i_nparts;
			ss << k;
			ss >> i_nparts;
			ss.clear();

			string name_f = fname + "_" + i_nparts + "_of_" + s_nparts + "subdomains.msh";
			fstream test_out;
			test_out.open(name_f.c_str(), ios::out | ios::trunc);

			Node *nod = 0;
			//GMSH test_out<<"$NOD"<<endl;
			//GMSH test_out<<(long)nodes_dom.size()<<endl;
			test_out << "#0#0#0#1#0.0#0#################################################################" << endl;
			test_out << "0 " << (long) nodes_dom.size() << " " << (long) eles_dom.size() << endl;
			for (j = 0; j < (long) nodes_dom.size(); j++)
			{
				nod = node_vector[nodes_dom[j]];
				//GMSH  test_out<<j+1<<"  "
				test_out << j << deli << nod->X() << deli << nod->Y() << deli << nod->Z() << endl;
			}
			//GMSH test_out<<"$ENDNOD"<<endl;
			//GMSH test_out<<"$ELE"<<endl;
			//GMSG test_out<<(long)eles_dom.size()<<endl;
			for (j = 0; j < (long) eles_dom.size(); j++)
			{
				ele = eles_dom[j];

				//GMSH  ele->WriteGmsh(test_out, k+1);

				test_out << j << deli << ele->getPatchIndex() << deli << ele->getName() << deli;
				for (kk = 0; kk < ele->getNodesNumber(useQuadratic); kk++)
					test_out << ele->getDomNodeIndex(kk) << deli;
				test_out << endl;
			}
			test_out.clear();
			test_out.close();
		}

	}
	part_out << "#STOP " << endl;
	part_out.clear();
	part_out.close();

	//
	delete ele_dom;
	//delete nod_dom;
}

void Mesh::findInternalNodes(const vector<long> &node_dom_idx, const int idom, vector<bool> &sdom_marked, vector<Node*> &internal_quad_nodes, vector<Node*> &internal_nodes)
{
	for (size_t j = 0; j < node_dom_idx.size(); j++)
	{
		if (node_dom_idx[j] == idom && !sdom_marked[j])
		{
			if (j >= (size_t)NodesNumber_Linear)
			{
				internal_quad_nodes.push_back(node_vector[j]);
			}
			internal_nodes.push_back(node_vector[j]);
			sdom_marked[j] = true; // avoid other subdomain use this node
		}
	}
}

int Mesh::readMaterialDataFile(const std::string &fpath, const std::string &mat_fname, vector<string> &m_datanames, int num_data, vector<string> &m_headers, std::vector<size_t> &m_header_marker_per_data, std::vector<double> &m_ele_val)
{
	string line_buffer;
	string mat_fname_abs = fpath + mat_fname;
	ifstream is_mat(mat_fname_abs.c_str());
	if (!is_mat.good())
	{
		cout << "Material data file " << mat_fname_abs << " does not exist" << endl;
		exit(1);
	}
	is_mat >> num_data;
	m_datanames.resize(num_data);
	m_header_marker_per_data.resize(num_data + 1);
	m_header_marker_per_data[0] = 0;
	for (int k = 0; k < num_data; k++)
	{
		string data_name;
		is_mat >> m_datanames[k];
		m_datanames[k] = fpath + m_datanames[k];
	}
	is_mat.close();
	// Read each data file
	for (int k = 0; k < num_data; k++)
	{
		is_mat.open(m_datanames[k].c_str());
		if (!is_mat.good())
		{
			cout << "Material data file " << m_datanames[k] << " does not exist" << endl;
			exit(1);
		}
		while (!is_mat.eof())
		{
			getline(is_mat, line_buffer);
			if (line_buffer.find("$DATA") != string::npos)
			{
				m_headers.push_back(line_buffer);
				const size_t ne = elem_vector.size();
				for (size_t ie = 0; ie < ne; ie++)
				{
					long index;
					double m_val;
					is_mat >> index >> m_val;
					//m_ele_idx.push_back(index);
					m_ele_val.push_back(m_val);
				}
			}
			else if (line_buffer.find("#STOP") != string::npos)
			{
				break;
			}
			else if (line_buffer.size() > 0)
			{
				m_headers.push_back(line_buffer);
			}
		}
		m_header_marker_per_data[k + 1] = m_headers.size();
		is_mat.clear();
		is_mat.close();
	}
	return num_data;
}


void Mesh::outputRenumedVTK(const std::string fname, const std::string s_nparts, const int num_parts, std::vector<long> &nnodes_sdom_start, long end, std::vector<Node*> &sbd_nodes, const bool is_quad, const bool osdom)
{
	std::string f_iparts = fname + "_renum_" + s_nparts + ".msh";
	ofstream os(f_iparts.c_str(), ios::out | ios::trunc);
	//os.setf(ios::scientific, std::ios::floatfield);
	os.precision(15);
	// Output renumbered mesh
	os << "#FEM_MSH\n";
	os << " $PCS_TYPE\n";
	os << "  NO_PCS\n";
	os << " $NODES\n  " << NodesNumber_Quadratic << endl;
	std::vector<Node*> vec_sorted_nodes(NodesNumber_Quadratic);
	for (int idom = 0; idom < num_parts; idom++)
	{
		const long start = nnodes_sdom_start[idom];
		//      const long end = nnodes_sdom_linear_elements[idom];
		if (idom < num_parts - 1)
			end = nnodes_sdom_start[idom + 1];
		else
			end = static_cast<long>(sbd_nodes.size());

		for (long i = start; i < end; i++)
		{
			Node* a_node = sbd_nodes[i];
			a_node->index = a_node->local_index;
			vec_sorted_nodes[a_node->index] = a_node;
			//a_node->Write(os);
		}
	}
	for (size_t i = 0; i < vec_sorted_nodes.size(); i++)
		vec_sorted_nodes[i]->Write(os);
	sbd_nodes.clear();
	os << " $ELEMENTS\n  " << elem_vector.size() << endl;
	for (size_t e = 0; e < elem_vector.size(); e++)
	{
		elem_vector[e]->WriteGSmsh(os, is_quad);
	}
	os << "#STOP" << endl;
	os.close();
	if (osdom)
	{
		//-----------------------------------------------------------
		/// VTK output
		// Elements in this subdomain
		f_iparts = fname + "_renum_" + s_nparts + ".vtk";
		//f_iparts = fname+"_"+str_buf+".vtk";
		ofstream os(f_iparts.c_str(), ios::out | ios::trunc);
		WriteVTK_Nodes(os, vec_sorted_nodes, 0);
		WriteVTK_Elements_of_Subdomain(os, elem_vector, 0, 0);
		//-----------------------------------------------------------
	}
}

void Mesh::findElementsInSubDomain(const vector<Node*>& internal_nodes, vector<Elem*>& subdom_internal_elements, vector<Elem*>& subdom_ghost_elements)
{
	for (size_t j = 0; j < internal_nodes.size(); j++)
	{
		Node* a_node = internal_nodes[j];
		// Search the elements connected to this nodes
		const size_t ne_rel = a_node->ElementsRelated.size();
		for (size_t k = 0; k < ne_rel; k++)
		{
			Elem* a_elem = elem_vector[a_node->ElementsRelated[k]];
			if (a_elem->getStatus())
				continue;

			a_elem->Marking(true);
			vector<int> ng_nodes; // non ghost nodes in ghost elements
			vector<int> g_nodes; // ghost nodes in ghost elements
			vector<int> g_nodes_L; // ghost linear nodes in ghost elements
			vector<int> ngl_nodes; // ghost nodes in ghost elements
			for (int kk = 0; kk < a_elem->getNodesNumber(useQuadratic); kk++)
			{
				if (a_elem->getNode(kk)->getStatus()) {
					ng_nodes.push_back(kk);
					if (a_elem->getNode(kk)->isQuadratic())
						ngl_nodes.push_back(kk);
				} else {
					g_nodes.push_back(kk);
					if (a_elem->getNode(kk)->isQuadratic())
						g_nodes_L.push_back(kk);
				}
			}
			// All nodes of this element are inside this subdomain
			if (g_nodes.empty())
			{
				subdom_internal_elements.push_back(a_elem);
			}
			else if (g_nodes.size() != static_cast<size_t>(a_elem->getNodesNumber(useQuadratic)))
			{
				subdom_ghost_elements.push_back(a_elem);
				// set ghost nodes
				a_elem->nnodes_gl = g_nodes_L.size();
				a_elem->ghost_nodes.resize(g_nodes.size());
				for (int kk = 0; kk < g_nodes.size(); kk++)
					a_elem->ghost_nodes[kk] = g_nodes[kk];
//				const int nn_gl = static_cast<int>(ng_nodes.size());
//				a_elem->nnodes_gl = ngl_nodes.size();
//				a_elem->ghost_nodes.resize(nn_gl);
//				for (int kk = 0; kk < nn_gl; kk++)
//					a_elem->ghost_nodes[kk] = ng_nodes[kk];
			}
		}
	}
}

void Mesh::findGhostNodesInSubDomain(const vector<Elem*>& subdom_ghost_elements, const bool is_quad, vector<Node*>& dom_ghost_linear_nodes, vector<Node*>& dom_ghost_quad_nodes)
{
	for (size_t j = 0; j < subdom_ghost_elements.size(); j++)
	{
		Elem* a_elem = subdom_ghost_elements[j];
		for (int k = 0; k < a_elem->getNodesNumber(is_quad); k++)
			a_elem->nodes[k]->Marking(false);
		// mark ghost nodes
		for (size_t k = 0; k < a_elem->ghost_nodes.size(); k++)
			a_elem->nodes[a_elem->ghost_nodes[k]]->Marking(true); //ghost_nodes actually hold internal nodes
	}
	//
	for (size_t j = 0; j < subdom_ghost_elements.size(); j++)
	{
		Elem* a_elem = subdom_ghost_elements[j];
		for (int k = 0; k < a_elem->getNodesNumber(is_quad); k++)
		{
			Node* a_node = a_elem->nodes[k];
			if (!a_node->getStatus())
				// ghost nodes are unmarked
				continue;

			a_node->Marking(true);
			if (k < a_elem->getNodesNumber(false))
				dom_ghost_linear_nodes.push_back(a_node);
			else
				dom_ghost_quad_nodes.push_back(a_node);
		}
	}
}

long Mesh::addSubDomainNodes(long node_id_shift, const vector<Node*>& internal_nodes, const vector<Node*>& internal_quad_nodes, const vector<Node*>& dom_ghost_linear_nodes, const vector<Node*>& dom_ghost_quad_nodes, vector<Node*>& sbd_nodes)
{
	// make a list of domain nodes
	long new_node_idx = 0; //node_id_shift;
	// add internal linear
	for (size_t j = 0; j < internal_nodes.size() - internal_quad_nodes.size(); j++)
	{
		Node* a_node = internal_nodes[j];
		a_node->index = new_node_idx++; //local node id
		sbd_nodes.push_back(a_node);
	}
	// add ghost linear
	for (size_t j = 0; j < dom_ghost_linear_nodes.size(); j++)
	{
		Node* a_node = dom_ghost_linear_nodes[j];
		a_node->index = new_node_idx++;
		sbd_nodes.push_back(a_node);
	}
	// add internal quad
	for (size_t j = 0; j < internal_quad_nodes.size(); j++)
	{
		Node* a_node = internal_quad_nodes[j];
		a_node->index = new_node_idx++;
		sbd_nodes.push_back(a_node);
	}
	// add ghost quad
	for (size_t j = 0; j < dom_ghost_quad_nodes.size(); j++)
	{
		Node* a_node = dom_ghost_quad_nodes[j];
		a_node->index = new_node_idx++;
		sbd_nodes.push_back(a_node);
	}
	return new_node_idx;
}

void Mesh::writeMatData(int num_data, const vector<string>& m_datanames, const string& dom_str, const vector<size_t>& m_header_marker_per_data, const vector<string>& m_headers, const long n_all_elements, const long nei_size, const vector<Elem*>& subdom_internal_elements, const std::string& deli, const vector<double>& m_ele_val, const vector<Elem*>& subdom_ghost_elements)
{
	ofstream os_mat;
	for (int mm = 0; mm < num_data; mm++)
	{
		string mat_ofile_name = m_datanames[mm] + dom_str;
		os_mat.open(mat_ofile_name.c_str(), ios::trunc);
		for (size_t mh = m_header_marker_per_data[mm]; mh < m_header_marker_per_data[mm + 1]; mh++)
		{
			os_mat << m_headers[mh] << endl;
		}
		const long e_shift = n_all_elements * mm;
		for (long j = 0; j < nei_size; j++)
		{
			const long entry_index = subdom_internal_elements[j]->getIndex() + e_shift;
			os_mat << j << deli << m_ele_val[entry_index] << endl;
		}
		for (size_t j = 0; j < subdom_ghost_elements.size(); j++)
		{
			const long entry_index = subdom_ghost_elements[j]->getIndex() + e_shift;
			os_mat << j + nei_size << deli << m_ele_val[entry_index] << endl;
		}
		os_mat << "#STOP" << endl;
		os_mat.clear();
		os_mat.close();
	}
}

void Mesh::outputSubDomainVTK(int idom, const std::string &fname, const std::string &dom_str, const std::string &s_nparts,
		std::vector<Node*> &sbd_nodes, std::vector<Elem*> &subdom_internal_elements, std::vector<Elem*> &subdom_ghost_elements,
		long nnodes_previous_sdom, long node_id_shift,int num_data,
		std::vector<std::string> &m_headers, std::vector<size_t> &m_header_marker_per_data,
		long n_all_elements, long nei_size, std::vector<double> &m_ele_val)
{
	//-----------------------------------------------------------
	// VTK output
	// Elements in this subdomain
	std::string f_iparts = fname + "_" + dom_str + "_of_" + s_nparts + "_subdomains.vtk";
	//f_iparts = fname+"_"+str_buf+".vtk";
	ofstream os(f_iparts.c_str(), ios::out | ios::trunc);

	WriteVTK_Nodes(os, sbd_nodes, nnodes_previous_sdom);
	WriteVTK_Elements_of_Subdomain(os, subdom_internal_elements, idom + 1, node_id_shift);

	/// Material data partitioning
	if (num_data > 0)
	{
		for (int mm = 0; mm < num_data; mm++)
		{
			// Partition
			os << "SCALARS " << m_headers[m_header_marker_per_data[mm] + 4] << " double 1\nLOOKUP_TABLE default" << endl;
			const long e_shift = n_all_elements * mm;
			for (long i = 0; i < nei_size; i++)
			{
				const long entry_index = subdom_internal_elements[i]->getIndex() + e_shift;
				os << m_ele_val[entry_index] << endl;
			}
		}
	}
	os.clear();
	os.close();

	//// Ghost elements in this subdomain
	f_iparts = fname + "_" + dom_str + "_ghost_of_" + s_nparts + "_subdomains.vtk";
	////f_iparts = fname+"_"+str_buf+"ghost.vtk";
	os.open(f_iparts.c_str(), ios::out | ios::trunc);
	WriteVTK_Nodes(os, sbd_nodes, nnodes_previous_sdom);
	WriteVTK_Elements_of_Subdomain(os, subdom_ghost_elements, 0, node_id_shift);
	if (num_data > 0)
	{
		for (int mm = 0; mm < num_data; mm++)
		{
			// Partition
			os << "SCALARS " << m_headers[m_header_marker_per_data[mm] + 4] << " double 1\nLOOKUP_TABLE default" << endl;
			const long e_shift = n_all_elements * mm;
			for (size_t i = 0; i < subdom_ghost_elements.size(); i++)
			{
				const long entry_index = subdom_ghost_elements[i]->getIndex() + e_shift;
				os << m_ele_val[entry_index] << endl;
			}
		}
	}
	os.clear();
	os.close();
	//-----------------------------------------------------------
}

/*!
 \brief void Mesh::ConstructSubDomain_by_Nodes

 Partition a mesh ny nodes

 02.2012 WW
 */

void Mesh::ConstructSubDomain_by_Nodes(const string fname, const string fpath, const std::string mat_fname, const int num_parts, const bool is_quad, const bool outut_subdomains)
{
	const long n_all_nodes = static_cast<long>(node_vector.size());
	const long nnodes_in_usage = useQuadratic ? NodesNumber_Quadratic : NodesNumber_Linear;

	// -------------------------------------------------------------------------
	// Read material data
	// -------------------------------------------------------------------------
	int num_data = 0;
	vector<string> m_headers;
	vector<size_t> m_header_marker_per_data;
	vector<string> m_datanames;
	vector<double> m_ele_val;
	if (mat_fname.size() != 0)
		num_data = readMaterialDataFile(fpath, mat_fname, m_datanames, num_data, m_headers, m_header_marker_per_data, m_ele_val);


	// -------------------------------------------------------------------------
	// Read .npart file
	// -------------------------------------------------------------------------
	vector<bool> vec_node_dom_marked(n_all_nodes);
	vector<long> vec_node_dom_idx(nnodes_in_usage);
	if (num_parts == 1) {
		for(long i=0; i<n_all_nodes; i++)
		{
			vec_node_dom_idx[i] = 0;
			vec_node_dom_marked[i] = false;
		}
	} else {
		read_npart_file(num_parts, fname, vec_node_dom_idx, vec_node_dom_marked);
	}

	// -------------------------------------------------------------------------
	// open partitioned file and write header
	// -------------------------------------------------------------------------
	const string s_nparts(number2str(num_parts));
	const string name_f = fname + "_partitioned_" + s_nparts + ".msh";
	fstream os_subd(name_f.c_str(), ios::out | ios::trunc);
	const std::string str_header = "Subdomain mesh "
			"(Domain nodes(quad); Domain nodes(linear); Inner elements; Ghost elements; Internal nodes(linear); Internal nodes(quad)) "
			"Global nodes(linear); Global nodes(quad); Global elements; "
			"Total integer variables of elements;Total integer variables of ghost elements  ";
//   name_f = "Subdomain mesh "
//           "(Nodes;  Nodes_linear; Elements; Ghost elements; Nodes of Linear elements; Nodes of quadratic elements) "
//            "Nodes of Linear whole elements; Nodes of whole quadratic elements; "
//       "Total integer variables of elements;Total integer variables of ghost elements  ";
	os_subd << str_header << endl;
	os_subd << num_parts;
	if (this->axisymmetry)
		os_subd << " AXISYMMETRY";
	os_subd << endl;
	setw(14);
	os_subd.precision(14);
	//os_subd.setf(ios::fixed, ios::scientific);
	os_subd.setf(ios::scientific);

	// -------------------------------------------------------------------------
	// collect sub-domain nodes
	// -------------------------------------------------------------------------
	vector<long> nnodes_sdom_start(num_parts);
	vector<long> nnodes_sdom_linear_elements(num_parts);
	vector<long> nnodes_sdom_quadratic_elements(num_parts);
	vector<size_t> position_node_file(num_parts);
	const long n_all_elements = static_cast<long>(elem_vector.size());
	vector<Node*> sbd_nodes;
	vector<Node*> vec_linear_nodes;
	vector<Node*> vec_quad_nodes;

	for (long i=0; i<n_all_nodes; i++)
	{
		this->node_vector[i]->original_index = this->node_vector[i]->index;
		if (this->node_vector[i]->index < this->NodesNumber_Linear)
			this->node_vector[i]->isQuadratic(true);
	}

	long node_id_shift = 0;
	long nnodes_previous_sdom = 0;
	for (int idom = 0; idom < num_parts; idom++)
	{
		cout << "Process partition: " << idom << endl;
		nnodes_sdom_start[idom] = nnodes_previous_sdom;

		long size_sbd_nodes = 0;
		vector<Node*> internal_nodes; // both linear and quad
		vector<Node*> internal_quad_nodes; // only quad
		vector<Node*> dom_ghost_linear_nodes, dom_ghost_quad_nodes;
		vector<Elem*> subdom_internal_elements;
		vector<Elem*> subdom_ghost_elements;
		if (num_parts==1) {
			for (long j = 0; j < n_all_nodes; j++) {
				internal_nodes.push_back(node_vector[j]);
				if (node_vector[j]->index >= NodesNumber_Linear)
					internal_quad_nodes.push_back(node_vector[j]);
				sbd_nodes.push_back(node_vector[j]);
				internal_nodes[j]->local_index = j;
				internal_nodes[j]->subdom_id = idom;
			}
			for (long j = 0; j < n_all_elements; j++)
				subdom_internal_elements.push_back(elem_vector[j]);
			size_sbd_nodes = n_all_nodes;
		} else {
			// Un-making all nodes and elements of the whole mesh
			for (long j = 0; j < n_all_nodes; j++)
				node_vector[j]->Marking(false);
			for (long j = 0; j < n_all_elements; j++)
				elem_vector[j]->Marking(false);

			// Find and mark internal nodes
			findInternalNodes(vec_node_dom_idx, idom, vec_node_dom_marked, internal_quad_nodes, internal_nodes);
			nnodes_sdom_linear_elements[idom] = static_cast<long>(internal_nodes.size() - internal_quad_nodes.size());
			nnodes_sdom_quadratic_elements[idom] = static_cast<long>(internal_nodes.size());
			for (size_t j = 0; j < internal_nodes.size(); j++)
			{
				Node* a_node = internal_nodes[j];
				a_node->Marking(true);
				a_node->local_index = j + node_id_shift; //internal node id should be continuous
				a_node->subdom_id = idom;

				if (internal_nodes[j]->isQuadratic())
					vec_linear_nodes.push_back(internal_nodes[j]);
				else
					vec_quad_nodes.push_back(internal_nodes[j]);
			}

			// Find elements in this domain
			findElementsInSubDomain(internal_nodes, subdom_internal_elements, subdom_ghost_elements);

			// Find ghost nodes
			findGhostNodesInSubDomain(subdom_ghost_elements, is_quad, dom_ghost_linear_nodes, dom_ghost_quad_nodes);

			// make a list of domain nodes
			addSubDomainNodes(node_id_shift, internal_nodes, internal_quad_nodes, dom_ghost_linear_nodes, dom_ghost_quad_nodes, sbd_nodes);
			size_sbd_nodes = static_cast<long>(sbd_nodes.size()) - nnodes_previous_sdom;

		}

		long nmb_element_idxs = 0, nmb_element_idxs_g = 0;
		// Count the total integer variables of this subdomain
		nmb_element_idxs = 4 * subdom_internal_elements.size(); // global id, mat id, ele type, nr. nodes
		for (size_t j = 0; j < subdom_internal_elements.size(); j++)
			nmb_element_idxs += subdom_internal_elements[j]->getNodesNumber(is_quad);
		nmb_element_idxs_g = 6 * subdom_ghost_elements.size(); // above + nr. linear ghost nodes, nr. total ghost nodes
		for (size_t j = 0; j < subdom_ghost_elements.size(); j++)
		{
			nmb_element_idxs_g += subdom_ghost_elements[j]->getNodesNumber(is_quad);
			nmb_element_idxs_g += static_cast<long>(subdom_ghost_elements[j]->ghost_nodes.size());
		}

		//----------------------------------------------------------------------------------
		// write into a file
		//----------------------------------------------------------------------------------
		string dom_str(number2str(idom));
		const std::string deli(" ");
		os_subd << size_sbd_nodes << deli << size_sbd_nodes - internal_quad_nodes.size() - dom_ghost_quad_nodes.size()
				<< deli << subdom_internal_elements.size() << deli << subdom_ghost_elements.size() << deli
				<< internal_nodes.size() - internal_quad_nodes.size() << deli << internal_nodes.size() << deli
				<< NodesNumber_Linear << deli << NodesNumber_Quadratic << deli << this->elem_vector.size() << deli
				<< nmb_element_idxs << deli << nmb_element_idxs_g << endl;

		position_node_file[idom] = os_subd.tellp();
		//os_subd<<"Nodes"<<endl;
		for (long j = 0; j < size_sbd_nodes; j++)
			sbd_nodes[j + nnodes_previous_sdom]->WriteWithEqsID(os_subd, is_quad);

		//os_subd<<"Elements"<<endl;
		const long nei_size = static_cast<long>(subdom_internal_elements.size());
		for (long j = 0; j < nei_size; j++)
			subdom_internal_elements[j]->WriteSubDOM(os_subd, node_id_shift, is_quad);

		//os_subd<<"Ghost elements"<<endl;
		for (size_t j = 0; j < subdom_ghost_elements.size(); j++)
		{
			Elem* a_elem = subdom_ghost_elements[j];
			a_elem->WriteSubDOM(os_subd, node_id_shift, is_quad);
			const int ngh_nodes = static_cast<int>(a_elem->ghost_nodes.size());

			os_subd << a_elem->nnodes_gl << deli << ngh_nodes << deli;
			for (int kk = 0; kk < ngh_nodes; kk++)
			{
				os_subd << a_elem->ghost_nodes[kk] << deli;
			}
			os_subd << endl;
		}

		//----------------------------------------------------------------------------------
		// Material data partitioning
		if (num_data > 0)
			writeMatData(num_data, m_datanames, dom_str, m_header_marker_per_data, m_headers, n_all_elements, nei_size, subdom_internal_elements, deli, m_ele_val, subdom_ghost_elements);
		//----------------------------------------------------------------------------------

		if (outut_subdomains)
			outputSubDomainVTK(idom, fname, dom_str, s_nparts, sbd_nodes, subdom_internal_elements, subdom_ghost_elements, nnodes_previous_sdom, node_id_shift, num_data, m_headers, m_header_marker_per_data, n_all_elements, nei_size, m_ele_val);

		node_id_shift += internal_nodes.size();
		nnodes_previous_sdom = static_cast<long>(sbd_nodes.size());
	}

	//----------------------------------------------------------------------------------
	// re-number eqs index for linear case
	//----------------------------------------------------------------------------------
	size_t eqs_index = 0;
	for (size_t i=0; i<vec_linear_nodes.size(); i++) {
		vec_linear_nodes[i]->eqs_index = eqs_index++;
	}
	for (size_t i=0; i<vec_quad_nodes.size(); i++) {
		vec_quad_nodes[i]->eqs_index = eqs_index++;
	}
	if (is_quad)
	{
		// re-number eqs index for quadatic case
		size_t eqs_index_Q = 0;
		for (int idom = 0; idom < num_parts; idom++)
		{
			const long start = nnodes_sdom_start[idom];
			long end = 0;
			if (idom < num_parts - 1)
				end = nnodes_sdom_start[idom + 1];
			else
				end = static_cast<long>(sbd_nodes.size());

			for (long i = start; i < end; i++)
			{
				Node* a_node = sbd_nodes[i];
				if (a_node->subdom_id == idom)
					a_node->eqs_index_Q = eqs_index_Q++;
			}
		}
	}

	//----------------------------------------------------------------------------------
	// Rewrite nodes with new node index
	//----------------------------------------------------------------------------------
	long end = 0;
	for (int idom = 0; idom < num_parts; idom++)
	{
		const long start = nnodes_sdom_start[idom];
		if (idom < num_parts - 1)
			end = nnodes_sdom_start[idom + 1];
		else
			end = static_cast<long>(sbd_nodes.size());

		os_subd.seekp(position_node_file[idom]);
		for (long i = start; i < end; i++)
		{
			Node* a_node = sbd_nodes[i];
			a_node->WriteWithEqsID(os_subd, is_quad);
			//a_node->Write(os_subd);
		}
	}

	os_subd.clear();
	os_subd.close();

	outputRenumedVTK(fname, s_nparts, num_parts, nnodes_sdom_start, end, sbd_nodes, is_quad, outut_subdomains);
}

// 02.2012. WW
void Mesh::WriteVTK_Nodes(std::ostream& os)
{
	size_t i;
	Node *a_node = NULL;

	os << "# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n" << endl;
	os << "DATASET UNSTRUCTURED_GRID" << endl;
	os << "POINTS " << node_vector.size() << " double" << endl;
	setw(14);
	os.precision(14);
	for (i = 0; i < node_vector.size(); i++)
	{
		a_node = node_vector[i];
		os << a_node->X() << " " << a_node->Y() << " " << a_node->Z() << endl;
	}

}

// 03.2012. WW
void Mesh::WriteVTK_Nodes(std::ostream& os, std::vector<Node*>& nod_vec, const size_t start)
{
	size_t i;
	Node *a_node = NULL;

	os << "# vtk DataFile Version 4.0\nGrid Partition by WW \nASCII\n" << endl;
	os << "DATASET UNSTRUCTURED_GRID" << endl;
	os << "POINTS " << nod_vec.size() - start << " double" << endl;
	setw(14);
	os.precision(14);
	for (i = start; i < nod_vec.size(); i++)
	{
		a_node = nod_vec[i];
		os << a_node->X() << " " << a_node->Y() << " " << a_node->Z() << endl;
	}

}

// 02.2012. WW
void Mesh::WriteVTK_Elements_of_Subdomain(std::ostream& os, std::vector<Elem*>& ele_vec, const int sbd_index, const long node_shift)
{
	size_t i;
	int j, k;
	int nne;

	j = 0;
	//-----------------------------------------------------------
	//  VTK output
	// Elements in this subdomain
	size_t ne0 = ele_vec.size();
	size_t size = ne0;

	string deli = " ";

	Elem *a_elem = NULL;

	for (i = 0; i < ne0; i++)
	{
		a_elem = ele_vec[i];
		nne = a_elem->getNodesNumber(useQuadratic);
		if (useQuadratic && a_elem->ele_Type == quadri)
			nne -= 1;

		size += nne;
	}
	os << "\nCELLS " << ne0 << deli << size << endl;

	// CELLs
	for (i = 0; i < ne0; i++)
	{
		a_elem = ele_vec[i];

		nne = a_elem->getNodesNumber(useQuadratic);
		if (useQuadratic && a_elem->ele_Type == quadri)
			nne -= 1;

		os << nne << deli;

		if (useQuadratic && a_elem->ele_Type == tet) // Tet
		{
			for (k = 0; k < 7; k++)
			{
				os << a_elem->nodes[k]->getIndex() - node_shift << deli;
			}

			for (k = 0; k < 3; k++)
			{
				j = (k + 2) % 3 + 7;
				os << a_elem->nodes[j]->getIndex() - node_shift << deli;
			}
		}
		else
		{
			for (k = 0; k < nne; k++)
			{
				os << a_elem->nodes[k]->getIndex() - node_shift << deli;
			}
		}

		os << endl;
	}
	os << endl;

	// CELL types
	os << "CELL_TYPES " << ne0 << endl;
	for (i = 0; i < ne0; i++)
	{
		a_elem = ele_vec[i];
		a_elem->WriteVTK_Type(os, useQuadratic);
	}
	os << endl;

	// Partition
	os << "CELL_DATA " << ne0 << endl;
	os << "SCALARS Partition int 1\nLOOKUP_TABLE default" << endl;
	for (i = 0; i < ne0; i++)
		os << sbd_index << endl;

}

void Mesh::WriteVTK_Elements(std::ostream& os)
{
	size_t ne0 = elem_vector.size();
	size_t size = ne0;

	const string deli = " ";

	for (auto a_elem : elem_vector)
	{
		auto nne = a_elem->getNodesNumber(useQuadratic);
		if (useQuadratic && a_elem->ele_Type == quadri)
			nne -= 1;

		size += nne;
	}

	os << "\nCELLS " << ne0 << deli << size << endl;

	// CELLs
	for (auto a_elem : elem_vector)
	{

		auto nne = a_elem->getNodesNumber(useQuadratic);
		if (useQuadratic && a_elem->ele_Type == quadri)
			nne -= 1;

		os << nne << deli;

		if (useQuadratic && a_elem->ele_Type == tet) // Tet
		{
			for (int k = 0; k < 7; k++)
				os << a_elem->nodes[k]->getIndex() << deli;

			for (int k = 0; k < 3; k++)
			{
				int j = 7 + k;
				//int j = (k + 2) % 3 + 7;
				//

				os << a_elem->nodes[j]->getIndex() << deli;
			}
		}
		else
		{
			for (int k = 0; k < nne; k++)
				os << a_elem->nodes[k]->getIndex() << deli;
		}

		os << endl;
	}
	os << endl;

	// CELL types
	os << "CELL_TYPES " << ne0 << endl;
	for (auto a_elem : elem_vector)
		a_elem->WriteVTK_Type(os, useQuadratic);

	os << endl;

//	// Partition
//	os << "CELL_DATA " << ne0 << endl;
//	os << "SCALARS Partition int 1\nLOOKUP_TABLE default" << endl;
//	for (size_t i = 0; i < ne0; i++)
//		os << sbd_index << endl;

	os << "CELL_DATA " << ne0 << endl;
	os << "SCALARS MatID int 1\nLOOKUP_TABLE default" << endl;
	for (auto a_elem : elem_vector)
		os << a_elem->getPatchIndex() << endl;

}

void Mesh::WriteVTK_Vertex(std::ostream& os)
{
	//-----------------------------------------------------------
	//  VTK output
	WriteVTK_Nodes(os);

	os << "\nCELLS " << node_vector.size() << " " << 2 * node_vector.size() << endl;

	// CELLs
	for (size_t i = 0; i < node_vector.size(); i++)
	{
		os << "1 " << i << "\n";
	}
	os << endl;

	// CELL types
	os << "CELL_TYPES " << node_vector.size() << endl;
	for (size_t i = 0; i < node_vector.size(); i++)
	{
		os << "1 \n";
	}
	os << endl;

	// Partition
	os << "POINT_DATA " << node_vector.size() << endl;
	os << "SCALARS Partition int 1\nLOOKUP_TABLE default" << endl;
	for (size_t i = 0; i < node_vector.size(); i++)
		os << "0\n";

}

void Mesh::Write2METIS(ostream& os)
{
	os << (long) elem_vector.size() << " ";

#ifdef METIS4_0
	int e_type =0;
	switch(elem_vector[0]->getElementType())
	{
		case line:
		cout<<"Not for 1D element"<<endl;
		exit(1);
		case quadri:
		e_type =4;
		break;
		case hex:
		e_type =3;
		break;
		case tri:
		e_type =1;
		break;
		case tet:
		e_type =2;
		break;
		case 6:
		cout<<"Not for prismal element"<<endl;
		abort();
	}
	os<<e_type;
#endif
	os << endl;
	for (long i = 0; i < (long) elem_vector.size(); i++)
	{
		elem_vector[i]->setOrder(useQuadratic);
		elem_vector[i]->Write_index(os);
	}
}

// Read grid for test purpose
void Mesh::ReadGrid(istream& is)
{
	long i, ne, nn, counter;
	int ibuff;
	double x, y, z;
	string buffer;
//   is.seekg(position);
	// Read description
	is >> buffer >> ws;
	// Read numbers of nodes and elements
	is >> ibuff >> nn >> ne >> ws;
	if (nn == 0 || ne == 0)
	{
		cout << "Error: number of elements or nodes is zero" << endl;
		exit(1);
	}

	// Read Nodes
	counter = 0;
	for (i = 0; i < nn; i++)
	{
		is >> ibuff >> x >> y >> z >> ws;
		Node* newNode = new Node(ibuff, x, y, z);
		newNode->Marking(true);
		node_vector.push_back(newNode);
		counter++;
	}
	if (counter != nn)
	{
		cout << "Error: number nodes do not match" << endl;
		exit(1);
	}
	NodesNumber_Linear = nn;
	NodesNumber_Quadratic = nn;

	// Read Elements
	counter = 0;
	for (i = 0; i < ne; i++)
	{
		Elem* newElem = new Elem(i);
		newElem->Read(is, this->node_vector, 1);
		newElem->Marking(true);
		elem_vector.push_back(newElem);
		counter++;
	}
	if (counter != ne)
	{
		cout << "Error: number elements do not match" << endl;
		exit(1);
	}

//   position = is.tellg();
}

void Mesh::ReadGridGeoSys(istream& is)
{
	string sub_line;
	string line_string;
	bool new_keyword = false;
	string hash("#");
	string sub_string, sub_string1;
	long i, ibuff;
	long no_elements;
	long no_nodes;
	double x, y, z;
	Node* newNode = NULL;
	Elem* newElem = NULL;
	//========================================================================
	// Keyword loop
	while (!new_keyword)
	{
		//if(!GetLineFromFile(line,fem_file))
		//  break;
		//line_string = line;
		getline(is, line_string);
		if (is.fail())
			break;
		/*
		 if(line_string.find(hash)!=string::npos)
		 {
		 new_keyword = true;
		 break;
		 }
		 */
		//....................................................................
		//....................................................................
		if (line_string.find("$AXISYMMETRY") != string::npos)   // subkeyword found
		{
			axisymmetry = true;
			continue;
		}
		//....................................................................
		if (line_string.find("$NODES") != string::npos)   // subkeyword found
		{
			is >> no_nodes >> ws;
			for (i = 0; i < no_nodes; i++)
			{
				is >> ibuff >> x >> y >> z >> ws;
				newNode = new Node(ibuff, x, y, z);
				node_vector.push_back(newNode);
			}
			continue;
		}
		//....................................................................
		if (line_string.find("$ELEMENTS") != string::npos)   // subkeyword found
		{
			is >> no_elements >> ws;
			for (i = 0; i < no_elements; i++)
			{
				newElem = new Elem(i);
				newElem->Read(is, this->node_vector, 0);
				newElem->Marking(true);
				elem_vector.push_back(newElem);
			}
			continue;
		}
	}
	//========================================================================
}

}   //end namespace

