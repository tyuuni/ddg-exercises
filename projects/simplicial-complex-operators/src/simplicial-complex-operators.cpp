// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    SparseMatrix<size_t> mat(mesh->nEdges(), mesh->nVertices());
    for (Edge e : mesh->edges()) {
        mat.insert(e.getIndex(), e.firstVertex().getIndex()) = 1;
        mat.insert(e.getIndex(), e.secondVertex().getIndex()) = 1;
    }
    return mat;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    SparseMatrix<size_t> mat(mesh->nFaces(), mesh->nEdges());
    for (Face f : mesh->faces()) {
        for (Edge e : f.adjacentEdges()) {
            mat.insert(f.getIndex(), e.getIndex()) = 1;
        }
    }
    return mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> vec(mesh->nVertices());
    for (size_t i = 0; i < mesh->nVertices(); i++) {
        vec[i] = subset.vertices.find(i) != subset.vertices.end();
    }
    return vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    Vector<size_t> vec(mesh->nEdges());
    for (size_t i = 0; i < mesh->nEdges(); i++) {
        vec[i] = subset.edges.find(i) != subset.edges.end();
    }
    return vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    Vector<size_t> vec(mesh->nFaces());
    for (size_t i = 0; i < mesh->nFaces(); i++) {
        vec[i] = subset.faces.find(i) != subset.faces.end();
    }
    return vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    MeshSubset ret(subset);
    for (auto it = subset.vertices.begin(); it != subset.vertices.end(); it++) {
        for (Edge edge : mesh->vertex(*it).adjacentEdges()) {
            ret.edges.insert(edge.getIndex());
        }
        for (Face face : mesh->vertex(*it).adjacentFaces()) {
            ret.faces.insert(face.getIndex());
        }
    }
    for (auto it = subset.edges.begin(); it != subset.edges.end(); it++) {
        Edge edge = mesh->edge(*it);
        for (Face face : edge.adjacentFaces()) {
            ret.faces.insert(face.getIndex());
        }
    }
    return ret;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    MeshSubset ret(subset);
    for (auto it = subset.edges.begin(); it != subset.edges.end(); it++) {
        Edge edge = mesh->edge(*it);
        ret.vertices.insert(edge.firstVertex().getIndex());
        ret.vertices.insert(edge.secondVertex().getIndex());
    }
    for (auto it = subset.faces.begin(); it != subset.faces.end(); it++) {
        Face face = mesh->face(*it);
        for (Edge edge : face.adjacentEdges()) {
            ret.edges.insert(edge.getIndex());
        }
        for (Vertex vertex : face.adjacentVertices()) {
            ret.vertices.insert(vertex.getIndex());
        }
    }
    return ret;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    MeshSubset s = star(closure(subset));
    MeshSubset ret = closure(star(subset));
    for (auto it = s.vertices.begin(); it != s.vertices.end(); it++) {
        ret.vertices.erase(*it);
    }
    for (auto it = s.edges.begin(); it != s.edges.end(); it++) {
        ret.edges.erase(*it);
    }
    for (auto it = s.faces.begin(); it != s.faces.end(); it++) {
        ret.faces.erase(*it);
    }
    return ret;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    for (auto it = subset.faces.begin(); it != subset.faces.end(); it++) {
        Face face = mesh->face(*it);
        for (Edge edge : face.adjacentEdges()) {
            if (subset.edges.find(edge.getIndex()) == subset.edges.end()) {
                return false;
            }
        }
    }
    for (auto it = subset.edges.begin(); it != subset.edges.end(); it++) {
        Edge edge = mesh->edge(*it);
        if (subset.vertices.find(edge.firstVertex().getIndex()) == subset.vertices.end() ||
            subset.vertices.find(edge.secondVertex().getIndex()) == subset.vertices.end()) {
            return false;
        }
    }
    return true;
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    if (!subset.faces.empty()) {
        if (isComplex(subset)) {
            for (auto it = subset.edges.begin(); it != subset.edges.end(); it++) {
                Edge edge = mesh->edge(*it);
                bool is_contained_in_subset = false;
                for (Face face : edge.adjacentFaces()) {
                    if (subset.faces.find(face.getIndex()) != subset.faces.end()) {
                        is_contained_in_subset = true;
                        break;
                    }
                }
                if (!is_contained_in_subset) {
                    return -1;
                }
            }
            return 2;
        }
        return -1;
    }
    if (!subset.edges.empty()) {
        if (isComplex(subset)) {
            return 1;
        }
        return -1;
    }
    if (!subset.vertices.empty()) {
        return 0;
    }
    return -1;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    MeshSubset ret;
    if (!subset.faces.empty()) {
        for (auto it = subset.edges.begin(); it != subset.edges.end(); it++) {
            Edge edge = mesh->edge(*it);
            int contained_count = 0;
            for (Face face : edge.adjacentFaces()) {
                if (subset.faces.find(face.getIndex()) != subset.faces.end()) {
                    contained_count++;
                    if (contained_count > 1) {
                        break;
                    }
                }
            }
            if (contained_count == 1) {
                ret.edges.insert(*it);
                ret.vertices.insert(edge.firstVertex().getIndex());
                ret.vertices.insert(edge.secondVertex().getIndex());
            }
        }
        return ret;
    }
    if (!subset.edges.empty()) {
        for (auto it = subset.vertices.begin(); it != subset.vertices.end(); it++) {
            Vertex vertex = mesh->vertex(*it);
            int contained_count = 0;
            for (Edge edge : vertex.adjacentEdges()) {
                if (subset.edges.find(edge.getIndex()) != subset.edges.end()) {
                    contained_count++;
                    if (contained_count > 1) {
                        break;
                    }
                }
            }
            if (contained_count == 1) {
                ret.vertices.insert(*it);
            }
        }
    }
    return ret;
}