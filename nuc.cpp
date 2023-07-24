#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <algorithm>
#include <vector>

struct Vertex
{
	std::vector<int> adj_cell_;
	std::vector<int> sub_ver_;
};

std::vector<int> compute_simple_bd(std::vector<int> mesh_tri, unsigned int facet_num)
{
	unsigned int n = *std::max_element(mesh_tri.begin(), mesh_tri.end())+1;
	std::vector<int> row_temp;
	row_temp.resize(n, 0);
	std::vector<std::vector<int> > amd;
	amd.resize(n, row_temp);

	// We collect the pair of points that need checking connectivity
	std::vector<std::pair<int, int> > P;
	P.clear();
	for (unsigned int i = 0; i < facet_num; ++i)
	{
		P.emplace_back(std::pair<int, int>(mesh_tri[i * 3 + 0], mesh_tri[i * 3 + 1]));
		P.emplace_back(std::pair<int, int>(mesh_tri[i * 3 + 1], mesh_tri[i * 3 + 2]));
		P.emplace_back(std::pair<int, int>(mesh_tri[i * 3 + 2], mesh_tri[i * 3 + 0]));
	}

	for (auto iter = P.begin(); iter != P.end(); ++iter)
	{
		int p1 = iter->first;
		int p2 = iter->second;
		if (amd[p2][p1] == 1)
			amd[p2][p1] = 0;
		else
			amd[p1][p2] = 1;
	}

	// Create the boundaries
	int count = 0;
	for (auto iter = amd.begin(); iter != amd.end(); ++iter)
	{
		for (auto iter2 = iter->begin(); iter2 != iter->end(); ++iter2)
		{
			if (*iter2 != 0)
				count++;
		}
	}

	std::vector<int> bd;
	// We find the first non-zero element
	int x = -1, y = -1;
	for (unsigned int i = 0; i < n; ++i)
	{
		for (unsigned int j = 0; j < n; ++j)
		{
			if (amd[i][j] != 0)
			{
				x = i;
				y = j;
				break;
			}
		}
		if (x != -1)
			break;
	}

	while (y != -1)
	{
		bd.emplace_back(x);
		amd[x][y] = 0;
		x = y;
		y = -1;
		for (unsigned int j = 0; j < n; ++j)
		{
			if (amd[x][j] != 0)
			{
				y = j;
				break;
			}
		}
	}

	return bd;
}


std::vector<std::vector<int> > compute_adjacency_matrix(std::vector<int> mesh_tri, unsigned int facet_num)
{
	std::vector<std::vector<int> > am;
	int max_num = *std::max_element(mesh_tri.begin(), mesh_tri.end())+1;
	std::vector<int> temp;
	temp.resize(max_num, 0);
	am.resize(max_num, temp);
	for (unsigned int i = 0; i < facet_num; ++i)
	{
		am[mesh_tri[i * 3 + 0]][mesh_tri[i * 3 + 1]] = 1;
		am[mesh_tri[i * 3 + 1]][mesh_tri[i * 3 + 0]] = 1;
		am[mesh_tri[i * 3 + 1]][mesh_tri[i * 3 + 2]] = 1;
		am[mesh_tri[i * 3 + 2]][mesh_tri[i * 3 + 1]] = 1;
		am[mesh_tri[i * 3 + 2]][mesh_tri[i * 3 + 0]] = 1;
		am[mesh_tri[i * 3 + 0]][mesh_tri[i * 3 + 2]] = 1;
	}

	return am;
}

void findAdjacentFacet(std::vector<int> mesh_tri, unsigned int facet_num, std::vector<std::vector<int> >& adjFacetList)
{
	int N = *std::max_element(mesh_tri.begin(), mesh_tri.end())+1;

	std::vector<int> row_temp;
	row_temp.resize(N, -1);
	std::vector<std::vector<int> > adj_matrix;
	adj_matrix.resize(N, row_temp);

	for (unsigned int i = 0; i < facet_num; ++i)
	{
		int c1 = mesh_tri[i * 3 + 0];
		int c2 = mesh_tri[i * 3 + 1];
		adj_matrix[c1][c2] = i;

		c1 = mesh_tri[i * 3 + 1];
		c2 = mesh_tri[i * 3 + 2];
		adj_matrix[c1][c2] = i;

		c1 = mesh_tri[i * 3 + 2];
		c2 = mesh_tri[i * 3 + 0];
		adj_matrix[c1][c2] = i;

	}

	std::vector<int> temp;
	adjFacetList.resize(facet_num, temp);
	for (unsigned int i = 0; i < facet_num; ++i)
	{
		int c1 = mesh_tri[i * 3 + 0];
		int c2 = mesh_tri[i * 3 + 1];
		adjFacetList[i].emplace_back(adj_matrix[c2][c1]);
		
		c1 = mesh_tri[i * 3 + 1];
		c2 = mesh_tri[i * 3 + 2];
		adjFacetList[i].emplace_back(adj_matrix[c2][c1]);

		c1 = mesh_tri[i * 3 + 2];
		c2 = mesh_tri[i * 3 + 0];
		adjFacetList[i].emplace_back(adj_matrix[c2][c1]);

	}

}

void edgeDivision(std::vector<double> mesh_ver, 
	std::vector<int> mesh_tri, unsigned int facet_num, 
	std::vector<double>& sub_division_ver, 
	std::vector<int>& sub_division_tri, 
	std::vector<Vertex>& vertices)
{
	vertices.clear();
	Vertex temp;
	vertices.resize(mesh_tri.size() / 3, temp);

	// Using the barocenter and the midpoint of edges, split original triangles

	std::vector<int> global_midpoint_index;
	std::vector<std::vector<int> > am = compute_adjacency_matrix(mesh_tri, facet_num);

	int n = 0; 
	for (unsigned int i = 0; i < am.size(); ++i)
	{
		for (unsigned int j = 0; j < am[0].size(); ++j)
		{
			if (am[i][j] != 0)
				n++;
		}
	}
	n /= 2;

	sub_division_ver.reserve(mesh_ver.size() + n * 3);
	sub_division_ver.assign(mesh_ver.begin(), mesh_ver.end());

	n = mesh_ver.size() / 3;
	for (unsigned int p = 0; p < facet_num; ++p)
	{
		// p is the index of triangles in the original mesh
		int i = mesh_tri[p * 3 + 0];
		int j = mesh_tri[p * 3 + 1];
		if (am[i][j] < 3)
		{
			double x = (mesh_ver[i * 3 + 0] + mesh_ver[j * 3 + 0]) / 2;
			double y = (mesh_ver[i * 3 + 1] + mesh_ver[j * 3 + 1]) / 2;
			double z = (mesh_ver[i * 3 + 2] + mesh_ver[j * 3 + 2]) / 2;
			sub_division_ver.emplace_back(x);
			sub_division_ver.emplace_back(y);
			sub_division_ver.emplace_back(z);
			am[i][j] = n;
			am[j][i] = n;
			n++;
			global_midpoint_index.emplace_back(n);
		}

		i = mesh_tri[p * 3 + 1];
		j = mesh_tri[p * 3 + 2];
		if (am[i][j] < 3)
		{
			double x = (mesh_ver[i * 3 + 0] + mesh_ver[j * 3 + 0]) / 2;
			double y = (mesh_ver[i * 3 + 1] + mesh_ver[j * 3 + 1]) / 2;
			double z = (mesh_ver[i * 3 + 2] + mesh_ver[j * 3 + 2]) / 2;
			sub_division_ver.emplace_back(x);
			sub_division_ver.emplace_back(y);
			sub_division_ver.emplace_back(z);
			am[i][j] = n;
			am[j][i] = n;
			n++;
			global_midpoint_index.emplace_back(n);
		}

		i = mesh_tri[p * 3 + 2];
		j = mesh_tri[p * 3 + 0];
		if (am[i][j] < 3)
		{
			double x = (mesh_ver[i * 3 + 0] + mesh_ver[j * 3 + 0]) / 2;
			double y = (mesh_ver[i * 3 + 1] + mesh_ver[j * 3 + 1]) / 2;
			double z = (mesh_ver[i * 3 + 2] + mesh_ver[j * 3 + 2]) / 2;
			sub_division_ver.emplace_back(x);
			sub_division_ver.emplace_back(y);
			sub_division_ver.emplace_back(z);
			am[i][j] = n;
			am[j][i] = n;
			n++;
			global_midpoint_index.emplace_back(n);
		}
	}

	// "am" now stores the indices of midpoint of all edges
	// Create the refined triangular mesh (for divisions)
	sub_division_tri.resize(mesh_tri.size()*4, 0);

	for (unsigned int i = 0; i < facet_num; ++i)
	{
		int verindex[4];
		verindex[0] = mesh_tri[i * 3 + 0];
		verindex[1] = mesh_tri[i * 3 + 1];
		verindex[2] = mesh_tri[i * 3 + 2];
		verindex[3] = mesh_tri[i * 3 + 0];

		// we calculate the barocenter
		double P[3][3];
		for (unsigned int j = 0; j < 3; ++j)
		{
			P[j][0] = mesh_ver[verindex[j] * 3 + 0];
			P[j][1] = mesh_ver[verindex[j] * 3 + 1];
			P[j][2] = mesh_ver[verindex[j] * 3 + 2];
		}

		double b[3];
		for(unsigned int j = 0; j < 3; ++j)
			b[j] = (P[0][j]+P[1][j]+P[2][j]) / 3;

		int b1 = sub_division_ver.size() / 3;
		sub_division_ver.emplace_back(b[0]);
		sub_division_ver.emplace_back(b[1]);
		sub_division_ver.emplace_back(b[2]);

		int v1 = verindex[0];
		int v2 = verindex[1];
		int v3 = verindex[2];
		int m1 = am[v1][v2];
		int m2 = am[v2][v3];
		int m3 = am[v3][v1];
		sub_division_tri[(i * 3 + 0) * 4 + 0] = b1;
		sub_division_tri[(i * 3 + 0) * 4 + 1] = m1;
		sub_division_tri[(i * 3 + 0) * 4 + 2] = v2;
		sub_division_tri[(i * 3 + 0) * 4 + 3] = m2;

		sub_division_tri[(i * 3 + 1) * 4 + 0] = b1;
		sub_division_tri[(i * 3 + 1) * 4 + 1] = m2;
		sub_division_tri[(i * 3 + 1) * 4 + 2] = v3;
		sub_division_tri[(i * 3 + 1) * 4 + 3] = m3;

		sub_division_tri[(i * 3 + 2) * 4 + 0] = b1;
		sub_division_tri[(i * 3 + 2) * 4 + 1] = m3;
		sub_division_tri[(i * 3 + 2) * 4 + 2] = v1;
		sub_division_tri[(i * 3 + 2) * 4 + 3] = m1;

		vertices[i].sub_ver_.emplace_back(i * 3 + 0);
		vertices[i].sub_ver_.emplace_back(i * 3 + 1);
		vertices[i].sub_ver_.emplace_back(i * 3 + 2);
	}
}


std::pair<int, int> locateFacet(std::vector<int> mesh_tri, int v1, int v2)
{
	std::pair<int, int>  result(-1, -1);
	
	for (unsigned int i = 0; i < mesh_tri.size() / 3; ++i)
	{
		if (mesh_tri[i * 3 + 0] == v1 && mesh_tri[i * 3 + 1] == v2)
		{
			result.first = i;
			result.second = 0;
			return result;
		}
		else if (mesh_tri[i * 3 + 1] == v1 && mesh_tri[i * 3 + 2] == v2)
		{
			result.first = i;
			result.second = 1;
			return result;
		}
		else if (mesh_tri[i * 3 + 2] == v1 && mesh_tri[i * 3 + 0] == v2)
		{
			result.first = i;
			result.second = 2;
			return result;
		}
	}
	return result;
}

std::vector<int> NUC(std::vector<int> mesh_tri, std::vector<double> mesh_ver, std::vector<Vertex> vertices)
{
	std::vector<int> result_path;

	std::vector<int> bd = compute_simple_bd(mesh_tri, mesh_tri.size()/3);
	if(bd.size() == 0)
	{
		std::cout << "The surface does not have a boundary. We return with empty result path. " << std::endl; 
		result_path.clear();
		return result_path;
	}

	std::vector<int> coverable_facets;
	coverable_facets.resize(mesh_tri.size() / 3, 1);

	// We find the first facet to be covered
	std::pair<int, int> init_temp = locateFacet(mesh_tri, bd[0], bd[1]);

	int init_facet_index = init_temp.first;

	if (init_temp.second == 0)
	{
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 0]);
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 1]);
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 2]);
	}
	else if (init_temp.second == 1)
	{
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 2]);
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 0]);
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 1]);
	}
	else if (init_temp.second == 2)
	{
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 1]);
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 2]);
		result_path.emplace_back(mesh_tri[init_facet_index * 3 + 0]);
	}
	
	std::vector<int> painted;
	painted.resize(vertices.size() * 3, 0);

	for (auto iter = result_path.begin(); iter != result_path.end(); ++iter)
	{
		painted[*iter] = 1;
	}

	// We assign connectivity of sub-facets
	std::vector<std::vector<int> > sub_adjacency_matrix;
	std::vector<int> temp;
	temp.resize(mesh_tri.size(), 0);
	sub_adjacency_matrix.resize( mesh_tri.size(), temp );

	unsigned int facet_num = mesh_tri.size() / 3;
	for (unsigned int i = 0; i < facet_num; ++i)
	{
		int v1 = mesh_tri[i * 3 + 0];
		int v2 = mesh_tri[i * 3 + 1];
		int v3 = mesh_tri[i * 3 + 2];

		// We deal with the first edge
		std::pair<int, int> temp = locateFacet(mesh_tri, v2, v1);
		int theOther = temp.first;
		if (theOther != -1 && coverable_facets[theOther] == 1)
		{
			int loc = temp.second;
			if (loc == 0)
			{
				sub_adjacency_matrix[3 * i + 2][3 * theOther + 0] = 1;
				sub_adjacency_matrix[3 * i + 0][3 * theOther + 2] = 1;
			}
			else if (loc == 1)
			{
				sub_adjacency_matrix[3 * i + 2][3 * theOther + 1] = 1;
				sub_adjacency_matrix[3 * i + 0][3 * theOther + 0] = 1;
			}
			else if (loc == 2)
			{
				sub_adjacency_matrix[3 * i + 2][3 * theOther + 2] = 1;
				sub_adjacency_matrix[3 * i + 0][3 * theOther + 1] = 1;
			}
		}

		// We deal with the second edge
		temp = locateFacet(mesh_tri, v3, v2);
		theOther = temp.first;
		if (theOther != -1 && coverable_facets[theOther] == 1)
		{
			int loc = temp.second;
			if (loc == 0)
			{
				sub_adjacency_matrix[3 * i + 0][3 * theOther + 0] = 1;
				sub_adjacency_matrix[3 * i + 1][3 * theOther + 2] = 1;
			}
			else if (loc == 1)
			{
				sub_adjacency_matrix[3 * i + 0][3 * theOther + 1] = 1;
				sub_adjacency_matrix[3 * i + 1][3 * theOther + 0] = 1;
			}
			else if (loc == 2)
			{
				sub_adjacency_matrix[3 * i + 0][3 * theOther + 2] = 1;
				sub_adjacency_matrix[3 * i + 1][3 * theOther + 1] = 1;
			}
		}

		// We deal with the third edge
		temp = locateFacet(mesh_tri, v1, v3);
		theOther = temp.first;
		if (theOther != -1 && coverable_facets[theOther] == 1)
		{
			int loc = temp.second;
			if (loc == 0)
			{
				sub_adjacency_matrix[3 * i + 1][3 * theOther + 0] = 1;
				sub_adjacency_matrix[3 * i + 2][3 * theOther + 2] = 1;
			}
			else if (loc == 1)
			{
				sub_adjacency_matrix[3 * i + 1][3 * theOther + 1] = 1;
				sub_adjacency_matrix[3 * i + 2][3 * theOther + 0] = 1;
			}
			else if (loc == 2)
			{
				sub_adjacency_matrix[3 * i + 1][3 * theOther + 2] = 1;
				sub_adjacency_matrix[3 * i + 2][3 * theOther + 1] = 1;
			}
		}
	}

	unsigned int i = 0;
	while (i < result_path.size())
	{
		auto& temp = sub_adjacency_matrix[result_path[i]];
		int loc = -1;
		for (auto iter = temp.begin(); iter != temp.end(); ++iter)
		{
			if (*iter == 1 && painted[iter-temp.begin()] == 0)
			{
				loc = iter - temp.begin();
				break;
			}
		}

		if (loc == -1)
		{
			i++;
			continue;
		}

		// By the index of loc we may identify the index of vertex
		int x = loc / 3;

		int adj_ver = x;

		std::vector<int> new_facets = vertices[adj_ver].sub_ver_;

		int b = std::find(new_facets.begin(), new_facets.end(), loc) - new_facets.begin();
		for (auto iter = new_facets.begin(); iter != new_facets.end(); ++iter)
			painted[*iter] = 1;

		std::vector<int> new_facets_temp;
		if (b != 0)
		{
			new_facets_temp.assign(new_facets.begin()+b, new_facets.end());
			new_facets_temp.insert(new_facets_temp.end(), new_facets.begin(), new_facets.begin() + b);
		}
		else
		{
			new_facets_temp.assign(new_facets.begin(), new_facets.end());
		}

		result_path.insert(result_path.begin() + i + 1, new_facets_temp.begin(), new_facets_temp.end());
		i++;
	}

	return result_path;
}

std::vector<double> evaluate_ours_in_saddle(std::vector<int> mesh_tri, std::vector<double> mesh_ver)
{
	// We do edge subdivision
	std::vector<std::vector<int> > adjFacetList;
	findAdjacentFacet(mesh_tri, mesh_tri.size() / 3, adjFacetList);

	std::vector<double> sub_division_ver;
	std::vector<int> sub_division_tri; // Each sub-facet has 4 vertices
	std::vector<Vertex> vertices;
	edgeDivision(mesh_ver, mesh_tri, mesh_tri.size() / 3, sub_division_ver, sub_division_tri, vertices);
	for (unsigned int i = 0; i < vertices.size(); ++i)
	{
		vertices[i].adj_cell_.assign(adjFacetList[i].begin(), adjFacetList[i].end());
	}

	std::vector<double> polygon_center; // center of each sub-facets
	// We find the center of sub-Facets
	polygon_center.resize((mesh_tri.size() / 3) * 3 * 3, 0);
	// One 3 is because each facet has been sub-divided into 3 sub-facets
	// Another 3 is because we are using a 1-dim vector to store a n*3 matrix
	for (unsigned int i = 0; i < (mesh_tri.size() / 3) * 3; ++i)
	{
		polygon_center[i * 3 + 0] = (sub_division_ver[sub_division_tri[i * 4 + 0] * 3 + 0] +
			sub_division_ver[sub_division_tri[i * 4 + 1] * 3 + 0] +
			sub_division_ver[sub_division_tri[i * 4 + 2] * 3 + 0] +
			sub_division_ver[sub_division_tri[i * 4 + 3] * 3 + 0]) / 4;
		polygon_center[i * 3 + 1] = (sub_division_ver[sub_division_tri[i * 4 + 0] * 3 + 1] +
			sub_division_ver[sub_division_tri[i * 4 + 1] * 3 + 1] +
			sub_division_ver[sub_division_tri[i * 4 + 2] * 3 + 1] +
			sub_division_ver[sub_division_tri[i * 4 + 3] * 3 + 1]) / 4;
		polygon_center[i * 3 + 2] = (sub_division_ver[sub_division_tri[i * 4 + 0] * 3 + 2] +
			sub_division_ver[sub_division_tri[i * 4 + 1] * 3 + 2] +
			sub_division_ver[sub_division_tri[i * 4 + 2] * 3 + 2] +
			sub_division_ver[sub_division_tri[i * 4 + 3] * 3 + 2]) / 4;
	}

	// We generate the template-independent motion
	std::vector<int> result_path = NUC(mesh_tri, mesh_ver, vertices);

	// We find the 3D coordinate of the skeleton
	std::vector<double> V;
	V.resize(result_path.size() * 3);
	for (unsigned int i = 0; i < result_path.size(); ++i)
	{
		V[i * 3 + 0] = polygon_center[result_path[i] * 3 + 0];
		V[i * 3 + 1] = polygon_center[result_path[i] * 3 + 1];
		V[i * 3 + 2] = polygon_center[result_path[i] * 3 + 2];
	}

	return V;
}


std::vector<double> add(std::vector<int> mesh_tri, std::vector<double> mesh_ver)
{
	std::vector<double> result;
	result = evaluate_ours_in_saddle(mesh_tri, mesh_ver);
    return result;
}

PYBIND11_MODULE(nuc, m)
{
    m.doc() = "pybind11 example";
    m.def("add", &add, "A function that adds two numbers.");
}


