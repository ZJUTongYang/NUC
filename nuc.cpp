#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <unordered_map>

namespace NUC
{
    void compute_adjacency_matrix_directed(const std::vector<int>& mesh_tri, int tri_num, int ver_num,
		std::unordered_map<int, int>& amd)
	{
		int index, first_ver, second_ver;
		for (unsigned int i = 0; i < tri_num; ++i)
		{
			if (mesh_tri[i * 3] == -1)
				continue;

			first_ver = mesh_tri[i * 3];
			second_ver = mesh_tri[i * 3 + 1];
			index = first_ver + ver_num * second_ver;
			amd[index] = i;

			first_ver = mesh_tri[i * 3 + 1];
			second_ver = mesh_tri[i * 3 + 2];
			index = first_ver + ver_num * second_ver;
			amd[index] = i;

			first_ver = mesh_tri[i * 3 + 2];
			second_ver = mesh_tri[i * 3];
			index = first_ver + ver_num * second_ver;
			amd[index] = i;
		}
	}

    struct Facet
    {
    public: 
		// The index of the facet in the mesh
        int index_;

        Facet(){}
        Facet(const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, 
            int tri_index, int subver_index, int child_in_parent_index)
        {
			child_.clear();
			index_ = tri_index;
			int v0 = mesh_tri[tri_index * 3];
			int v1 = mesh_tri[tri_index * 3 + 1];
			int v2 = mesh_tri[tri_index * 3 + 2];
			m_ = { (mesh_ver[v0 * 3] + mesh_ver[v1 * 3]) / 2, (mesh_ver[v0 * 3 + 1] + mesh_ver[v1 * 3 + 1]) / 2, (mesh_ver[v0 * 3 + 2] + mesh_ver[v1 * 3 + 2]) / 2 ,
				(mesh_ver[v1 * 3] + mesh_ver[v2 * 3]) / 2, (mesh_ver[v1 * 3 + 1] + mesh_ver[v2 * 3 + 1]) / 2, (mesh_ver[v1 * 3 + 2] + mesh_ver[v2 * 3 + 2]) / 2 ,
				(mesh_ver[v2 * 3] + mesh_ver[v0 * 3]) / 2, (mesh_ver[v2 * 3 + 1] + mesh_ver[v0 * 3 + 1]) / 2, (mesh_ver[v2 * 3 + 2] + mesh_ver[v0 * 3 + 2]) / 2
			};
			b_ = { (mesh_ver[v0 * 3] + mesh_ver[v1 * 3] + mesh_ver[v2 * 3]) / 3,
				(mesh_ver[v0 * 3 + 1] + mesh_ver[v1 * 3 + 1] + mesh_ver[v2 * 3 + 1]) / 3,
				(mesh_ver[v0 * 3 + 2] + mesh_ver[v1 * 3 + 2] + mesh_ver[v2 * 3 + 2]) / 3
			};
			polygon_center_ = { (mesh_ver[v1 * 3] + m_[1 * 3] + b_[0] + m_[0 * 3]) / 4, (mesh_ver[v1 * 3 + 1] + m_[1 * 3 + 1] + b_[0 + 1] + m_[0 * 3 + 1]) / 4, (mesh_ver[v1 * 3 + 2] + m_[1 * 3 + 2] + b_[0 + 2] + m_[0 * 3 + 2]) / 4,
								(mesh_ver[v2 * 3] + m_[2 * 3] + b_[0] + m_[1 * 3]) / 4, (mesh_ver[v2 * 3 + 1] + m_[2 * 3 + 1] + b_[0 + 1] + m_[1 * 3 + 1]) / 4, (mesh_ver[v2 * 3 + 2] + m_[2 * 3 + 2] + b_[0 + 2] + m_[1 * 3 + 2]) / 4,
								(mesh_ver[v0 * 3] + m_[0 * 3] + b_[0] + m_[2 * 3]) / 4, (mesh_ver[v0 * 3 + 1] + m_[0 * 3 + 1] + b_[0 + 1] + m_[2 * 3 + 1]) / 4, (mesh_ver[v0 * 3 + 2] + m_[0 * 3 + 2] + b_[0 + 2] + m_[2 * 3 + 2]) / 4
			};

			subver_index_ = subver_index;
			child_in_parent_index_ = child_in_parent_index;
        }
        void getCyclicCoveragePath(std::vector<int>& int_path, std::vector<double>& double_path) const
		{
			int subver_index = subver_index_ % 3;
			if (subver_index == 0)
			{
				int_path = { index_ * 3, index_ * 3 + 1, index_ * 3 + 2 };
				double_path.assign(polygon_center_.begin(), polygon_center_.end());
			}
			else if (subver_index == 1)
			{
				int_path = { index_ * 3 + 1, index_ * 3 + 2, index_ * 3 };
				double_path = { polygon_center_[3], polygon_center_[4], polygon_center_[5], polygon_center_[6], polygon_center_[7], polygon_center_[8], polygon_center_[0], polygon_center_[1], polygon_center_[2] };
			}
			else if (subver_index == 2)
			{
				int_path = { index_ * 3 + 2, index_ * 3, index_ * 3 + 1 };
				double_path = { polygon_center_[6], polygon_center_[7], polygon_center_[8], polygon_center_[0], polygon_center_[1], polygon_center_[2], polygon_center_[3], polygon_center_[4], polygon_center_[5] };
			}
		}

		std::vector<Facet*> child_;
		int child_in_parent_index_;
    private:

		std::vector<double> m_; // The Euclidean position of the midpoint of the three edges
		std::vector<double> b_; // The Euclidean position of the barocenter of the facet
		std::vector<double> polygon_center_; // The Euclidean position of the barocenter of sub-facets
		int subver_index_; // The order of the sub-facets designated by who the parent node is
    };

	Facet* insertTreeNode(Facet* parent, const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, 
		int tri_index, int subver_index, int child_in_parent_index)
	{
		parent->child_.emplace_back(new Facet(mesh_tri, mesh_ver, tri_index, subver_index, child_in_parent_index));
		return parent->child_.back();
	}

	void deleteTreeBranch(Facet* node)
	{
		if (node == nullptr)
			return;
		for (Facet* child : node->child_)
		{
			deleteTreeBranch(child);
		}
		delete node;
	}

	void traverse(Facet* node, const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, 
		std::vector<int>& topological_path, std::vector<double>& geometric_path, int position_to_insert)
	{
		if (node == nullptr)
			return;

		std::vector<int> cyclic_int_path;
		std::vector<double> cyclic_double_path;
		node->getCyclicCoveragePath(cyclic_int_path, cyclic_double_path);

		topological_path.insert(topological_path.begin() + position_to_insert, cyclic_int_path.begin(), cyclic_int_path.end());
		geometric_path.insert(geometric_path.begin() + position_to_insert * 3, cyclic_double_path.begin(), cyclic_double_path.end());

		int the_tri_index = node->index_;

		int v0 = mesh_tri[the_tri_index * 3];
		int v1 = mesh_tri[the_tri_index * 3+1];
		int v2 = mesh_tri[the_tri_index * 3+2];
		int child_in_parent_index;
		int loc;
		for (auto iter = node->child_.begin(); iter != node->child_.end(); ++iter)
		{
			child_in_parent_index = (*iter)->child_in_parent_index_;
			if (child_in_parent_index == 0)
			{
				loc = std::find(topological_path.begin(), topological_path.end(), the_tri_index * 3 + 2) - topological_path.begin();
			}
			else if (child_in_parent_index == 1)
			{
				loc = std::find(topological_path.begin(), topological_path.end(), the_tri_index * 3) - topological_path.begin();
			}
			else if (child_in_parent_index == 2)
			{
				loc = std::find(topological_path.begin(), topological_path.end(), the_tri_index * 3 + 1) - topological_path.begin();
			}

			traverse(*iter, mesh_tri, mesh_ver, topological_path, geometric_path, loc + 1);
		}
	}

	std::pair<std::vector<int>, std::vector<double> > nuc(const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver, int initial_tri_index)
    {
		std::vector<int> topological_coverage_path;
		std::vector<double> geometric_coverage_path;

		std::cout << "start generating non-revisiting uniform coverage path" << std::endl;
		Facet* root_ = nullptr;

		int tri_num = mesh_tri.size() / 3;
		int ver_num = mesh_ver.size() / 3;

		// We first collect the adjacency among facets
		std::unordered_map<int, int> amd; // <first_vertex_index + second_vertex_index * ver_num, facet_index>
		compute_adjacency_matrix_directed(mesh_tri, tri_num, ver_num, amd);

		std::cout << "size of amd: " << amd.size() << std::endl;
		// We avoid repetitive coverage
		std::vector<int> covered;
		covered.resize(tri_num, 0);

		// Invalid facets are marked as "covered", so that they are not considered during path deformation
		for (unsigned int i = 0; i < tri_num; ++i)
		{
			if (mesh_tri[i * 3] == -1)
			{
				covered[i] = 1;
			}
		}

		// The vertex indices of the source facet
		int v0 = mesh_tri[initial_tri_index * 3];
		int v1 = mesh_tri[initial_tri_index * 3 + 1];
		int v2 = mesh_tri[initial_tri_index * 3 + 2];

		// To avoid the problem in Theorem 1 in the paper, we select the best order of the sub-facets of the source facet
		// If any of the edges does not have adjacent facets, we select it as the "back"
		int subver_index;
		if (amd.find(v1 + v0 * ver_num) == amd.end())
		{
			subver_index = 0;
		}
		else if (amd.find(v2 + v1 * ver_num) == amd.end())
		{
			subver_index = 1;
		}
		else if (amd.find(v0 + v2 * ver_num) == amd.end())
		{
			subver_index = 2;
		}
		else
		{
			subver_index = 0;
		}
		
		// We put the initial facet at the root position
		root_ = new Facet(mesh_tri, mesh_ver, initial_tri_index, subver_index, -1);

		covered[initial_tri_index] = 1;
		std::list<Facet*> Q;
		Q.emplace_back(root_);

		int first_ver, second_ver;
		while (!Q.empty())
		{
			Facet* the_node = Q.front();
			Q.pop_front();

			int tri_index = the_node->index_;
			v0 = mesh_tri[tri_index * 3];
			v1 = mesh_tri[tri_index * 3 + 1];
			v2 = mesh_tri[tri_index * 3 + 2];

			first_ver = v0;
			second_ver = v1;
			std::unordered_map<int, int>::iterator iter = amd.find(second_ver + first_ver * ver_num);
			if (iter != amd.end())
			{
				int adj_tri_index = iter->second;
				if (covered[adj_tri_index] == 0)
				{
					if (mesh_tri[adj_tri_index * 3] == second_ver && mesh_tri[adj_tri_index * 3 + 1] == first_ver)
					{
						subver_index = 0;
					}
					else if (mesh_tri[adj_tri_index * 3 + 1] == second_ver && mesh_tri[adj_tri_index * 3 + 2] == first_ver)
					{
						subver_index = 1;
					}
					else if (mesh_tri[adj_tri_index * 3 + 2] == second_ver && mesh_tri[adj_tri_index * 3] == first_ver)
					{
						subver_index = 2;
					}
					else
					{
						subver_index = 0;
					}
					Facet* the_child = insertTreeNode(the_node, mesh_tri, mesh_ver, adj_tri_index, subver_index, 0);
					covered[adj_tri_index] = 1;
					Q.emplace_back(the_child);
				}
			}

			first_ver = v1;
			second_ver = v2;
			iter = amd.find(second_ver + first_ver * ver_num);
			if (iter != amd.end())
			{
				int adj_tri_index = iter->second;
				if (covered[adj_tri_index] == 0)
				{
					if (mesh_tri[adj_tri_index * 3] == second_ver && mesh_tri[adj_tri_index * 3 + 1] == first_ver)
					{
						subver_index = 0;
					}
					else if (mesh_tri[adj_tri_index * 3 + 1] == second_ver && mesh_tri[adj_tri_index * 3 + 2] == first_ver)
					{
						subver_index = 1;
					}
					else if (mesh_tri[adj_tri_index * 3 + 2] == second_ver && mesh_tri[adj_tri_index * 3] == first_ver)
					{
						subver_index = 2;
					}
					else
					{
						subver_index = 0;
					}
					Facet* the_child = insertTreeNode(the_node, mesh_tri, mesh_ver, adj_tri_index, subver_index, 1);
					covered[adj_tri_index] = 1;
					Q.emplace_back(the_child);
				}
			}

			first_ver = v2;
			second_ver = v0;
			iter = amd.find(second_ver + first_ver * ver_num);
			if (iter != amd.end())
			{
				int adj_tri_index = iter->second;
				if (covered[adj_tri_index] == 0)
				{
					if (mesh_tri[adj_tri_index * 3] == second_ver && mesh_tri[adj_tri_index * 3 + 1] == first_ver)
					{
						subver_index = 0;
					}
					else if (mesh_tri[adj_tri_index * 3 + 1] == second_ver && mesh_tri[adj_tri_index * 3 + 2] == first_ver)
					{
						subver_index = 1;
					}
					else if (mesh_tri[adj_tri_index * 3 + 2] == second_ver && mesh_tri[adj_tri_index * 3] == first_ver)
					{
						subver_index = 2;
					}
					else
					{
						subver_index = 0;
					}

					Facet* the_child = insertTreeNode(the_node, mesh_tri, mesh_ver, adj_tri_index, subver_index, 2);
					covered[adj_tri_index] = 1;
					Q.emplace_back(the_child);
				}
			}
		}

		// We report the geometric coverage path
		topological_coverage_path.clear();
		geometric_coverage_path.clear();

		int position_to_insert = 0;		
		traverse(root_, mesh_tri, mesh_ver, topological_coverage_path, geometric_coverage_path, position_to_insert);

		std::cout << "size of int path: " << topological_coverage_path.size() << std::endl;

		deleteTreeBranch(root_);

		return std::pair<std::vector<int>, std::vector<double> >(topological_coverage_path, geometric_coverage_path);
    }

    std::pair<std::vector<int>, std::vector<double> > nuc(const std::vector<int>& mesh_tri, const std::vector<double>& mesh_ver)
	{		
		// Here we allow for [-1, -1, -1] triangle facet, so we need to find the first valid facet
		int initial_tri_index = -1;
		int tri_num = mesh_tri.size()/3;
		
		for(unsigned int i = 0; i < tri_num; ++i)
		{
			if(mesh_tri[i*3] != -1)
			{
				initial_tri_index = i;
				break;
			}
		}

		if(initial_tri_index == -1)
			return std::pair<std::vector<int>, std::vector<double> >(std::vector<int>(), std::vector<double>());

		return nuc(mesh_tri, mesh_ver, initial_tri_index);
	}

}

PYBIND11_MODULE(nuc_tmech23, m)
{
	m.def("run", static_cast<std::pair<std::vector<int>, std::vector<double> > (*)(const std::vector<int>&, const std::vector<double>&)>(&NUC::nuc));
	m.def("run", static_cast<std::pair<std::vector<int>, std::vector<double> > (*)(const std::vector<int>&, const std::vector<double>&, int)>(&NUC::nuc));
}
