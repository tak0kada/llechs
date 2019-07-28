#pragma once

#include <array>
#include <vector>
#include <list>
#include <set>
#include <string>
#include <unordered_map>
#include <queue>
#include <utility>
#include <iostream>
#include <fstream>
#include <exception>
#include <boost/dynamic_bitset.hpp>
#include "happly/happly.h"

namespace std
{

template<>
struct hash<std::pair<std::size_t, std::size_t>>
{
    std::size_t operator()(const std::pair<std::size_t, std::size_t>& p) const
    {
        return std::hash<std::size_t>()(p.first) ^ std::hash<std::size_t>()(p.second);
    }
};

}


namespace post_acvd
{

struct Vertex
{
    std::size_t idx;
    double x;
    double y;
    double z;

    Vertex(const std::size_t& idx, const double& x, const double& y, const double& z)
    : idx{idx}, x{x}, y{y}, z{z}
    {}
};

struct Face
{
    std::size_t idx;
    std::array<std::size_t, 3> vertex;

    Face(const std::size_t& idx, const std::size_t& v0, const std::size_t& v1, const std::size_t& v2)
    : idx{idx}, vertex{v0, v1, v2}
    {}
};

struct DualGraph
{
private:
    std::unordered_map<std::size_t, std::list<std::size_t>> adj_list;

public:
    DualGraph(){};
    DualGraph(const std::size_t& nFace)
    : adj_list(nFace)
    {}

    void add_node(const std::size_t& i)
    {
        if (adj_list.find(i) == adj_list.end())
        {
            adj_list[i] = {};
        }
        else
        {
            std::cerr << "DualGraph already contains the node (face): " + std::to_string(i) + ".";
        }
    }
    void remove_node(const std::size_t& i)
    {
        if (auto cit = adj_list.find(i); cit == adj_list.end())
        {
            throw std::runtime_error("ERROR: DualGraph does not contain such node: " + std::to_string(i) + ".");
        }
        else
        {
            const auto& neighbor = cit->second;
            for (const auto& c: neighbor)
            {
                auto it = std::find(adj_list[c].begin(), adj_list[c].end(), i);
                adj_list[c].erase(it);
            }
            adj_list.erase(cit);
        }
    }

    void add_edge(const std::size_t& i, const std::size_t& j)
    {
        if (adj_list.find(i) == adj_list.end())
        {
            throw std::runtime_error("ERROR: DualGraph does not contain such node: " + std::to_string(i) + ".");
        }
        else
        {
            adj_list[i].push_back(j);
            adj_list[j].push_back(i);
        }
    }
    void remove_edge(const std::size_t& i, const std::size_t& j)
    {
        if (adj_list.find(i) == adj_list.end() || adj_list.find(j) == adj_list.end() ||
                std::find(adj_list[i].begin(), adj_list[i].end(), j) == adj_list[i].end())
        {
            throw std::runtime_error("ERROR: DualGraph does not contain such edge: " + std::to_string(i) + "-" + std::to_string(j) + ".");
        }
        else
        {
            auto it = std::find(adj_list[i].begin(), adj_list[i].end(), j);
            adj_list[i].erase(it);
            it = std::find(adj_list[j].begin(), adj_list[j].end(), i);
            adj_list[j].erase(it);
        }
    }

    std::vector<std::size_t> neighbor(const std::size_t& i) const
    {
        const auto it = adj_list.find(i);
        std::vector<std::size_t> ret(it->second.size());
        std::copy(it->second.begin(), it->second.end(), ret.begin());
        return ret;
    }
    int degree(const std::size_t& i) const
    {
        const auto it = adj_list.find(i);
        return it->second.size();
    }
    std::vector<std::size_t> allnode() const
    {
        std::vector<std::size_t> ret;
        ret.reserve(adj_list.size());
        for (const auto kv: adj_list)
        {
            ret.push_back(kv.first);
        }
        return ret;
    }
    std::size_t size() const
    {
        return adj_list.size();
    }
};

inline int degree(const std::size_t& i, const DualGraph& g)
{
    return g.degree(i);
}

// remove faces with edges that do not connect to any face
std::vector<Face> rm_border(const std::vector<Face> face)
{
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> conn(face.size() * 3 / 2 + 10);
    for (const auto& f: face)
    {
        std::array<std::size_t, 3> v = f.vertex;
        std::sort(v.begin(), v.end());
        conn[{v[0], v[1]}].push_back(f.idx);
        conn[{v[1], v[2]}].push_back(f.idx);
        conn[{v[0], v[2]}].push_back(f.idx);
    }

    boost::dynamic_bitset<> idx(face.size());
    while (true)
    {
        std::size_t n_rm{};
        for (auto& kv: conn)
        {
            if (kv.second.size() == 1)
            {
                kv.second.clear();
                idx[kv.second[0]] = true;
                ++n_rm;
            }
        }
        if (n_rm == 0)
        {
            break;
        }
    }

    std::vector<Face> ret;
    ret.reserve(face.size() - idx.count());
    for (std::size_t i = 0; i < face.size(); ++i)
    {
        if (idx[i] == false)
        {
            ret.push_back((face)[i]);
        }
    }
    return ret;
}

// void rm_iregular_node(DualGraph *g)
// {
//     std::vector<std::size_t> rm;
//     for (auto&& i: g->allnode())
//     {
//         if (g->degree(i) != 3)
//         {
//             rm.push_back(i);
//         }
//     }
//     for (auto&& i: rm)
//     {
//         g->remove_node(i);
//     }
// }

std::vector<std::size_t> find_debris(const DualGraph& g)
{
    const std::vector<std::size_t> node{g.allnode()};
    const std::size_t n_node = 1 + *std::max_element(node.begin(), node.end());
    boost::dynamic_bitset<> conn_default(n_node); //default: false
    boost::dynamic_bitset<> visited(n_node);
    conn_default.flip();
    visited.flip();
    for (const auto& n: node)
    {
        conn_default[n] = false;
        visited[n] = false;
    }

    boost::dynamic_bitset<> conn; //default: false
    while (conn.count() < g.size() * 0.9)
    {
        conn = conn_default;
        std::size_t start{};
        while (start < g.size())
        {
            if (visited[start] == false)
            {
                break;
            }
            ++start;
        }
        visited[start] = true;
        conn[start] = true;

        std::queue<std::size_t> q;
        q.push(start);
        while(!q.empty())
        {
            const std::size_t cur{q.front()};
            const std::vector<std::size_t> neighbor = g.neighbor(cur);
            if (neighbor.size() == 0)
            {
                break;
            }
            for (const std::size_t& n: neighbor)
            {
                if (!conn[n])
                {
                    q.push(n);
                    conn[n] = true;
                    visited[n] = true;
                }
            }
            q.pop();
        }
    }

    conn.flip();
    std::vector<std::size_t> ret;
    ret.reserve(conn.count());
    for (std::size_t pos = conn.find_first(); pos != conn.npos; pos = conn.find_next(pos))
    {
        ret.push_back(pos);
    }
    return ret;
}

void rm_debris(DualGraph *g)
{
    std::vector<std::size_t> rm;
    rm = find_debris(*g);
    for (const std::size_t& i: rm)
    {
        g->remove_node(i);
    }
}


DualGraph make_dual(const std::vector<Face>& face,
        const std::unordered_map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>>& conn)
{
    DualGraph graph(face.size());
    for (const Face& f: face)
    {
        graph.add_node(f.idx);
    }
    for (const auto& kv: conn)
    {
        const auto& v = kv.second;
        for (std::size_t i = 0; i < v.size() - 1; ++i)
        {
            for (std::size_t j = i + 1; j < v.size(); ++j)
            {
                graph.add_edge(v[i], v[j]);
            }
        }
    }
    return graph;
}

std::pair<std::vector<Vertex>, std::vector<Face>> parse_ply(const std::string& path)
{
    happly::PLYData plyIn(path);
    plyIn.validate();

    std::vector<std::array<double, 3>> v_tmp = plyIn.getVertexPositions();
    std::vector<std::vector<std::size_t>> f_tmp = plyIn.getFaceIndices();

    std::vector<Vertex> vertex;
    vertex.reserve(v_tmp.size());
    for (std::size_t vi = 0; vi < v_tmp.size(); ++vi)
    {
        vertex.emplace_back(vi, v_tmp[vi][0], v_tmp[vi][1], v_tmp[vi][2]);
    }
    std::vector<Face> face;
    face.reserve(f_tmp.size());
    for (std::size_t fi = 0; fi < f_tmp.size(); ++fi)
    {
        face.emplace_back(fi, f_tmp[fi][0], f_tmp[fi][1], f_tmp[fi][2]);
    }

    return {vertex, face};
}

std::pair<std::vector<Vertex>, std::vector<Face>> parse_obj(const std::string& path)
{
    std::ifstream ifs{path};
    if (ifs.fail())
    {
        throw std::runtime_error("ERROR: Cannot open obj file: " + path + ".");
    }

    std::vector<Vertex> vertex;
    std::vector<Face> face;
    std::string buf{};
    std::size_t vi{}, fi{};
    while (ifs >> buf)
    {
        if (buf == "v")
        {
            double x, y, z;
            ifs >> x >> y >> z;
            vertex.push_back({vi, x, y, z});
            ++vi;
        }
        else if (buf == "f")
        {
            std::size_t v0, v1, v2;
            ifs >> v0 >> v1 >> v2;
            face.push_back({fi, --v0, --v1, --v2});
            ++fi;
        }
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return {vertex, face};
}

std::pair<std::vector<std::array<double, 3>>, std::vector<std::array<std::size_t, 3>>> read_ply(const std::string& path)
{
    const auto [vertex_org, face_org] = parse_ply(path);
    const auto face = rm_border(face_org);

    // build edge-face coonection
    std::unordered_map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> conn(face.size() * 3 / 2 + 10);
    for (const auto& f: face)
    {
        std::array<std::size_t, 3> v = f.vertex;
        std::sort(v.begin(), v.end());
        conn[{v[0], v[1]}].push_back(f.idx);
        conn[{v[1], v[2]}].push_back(f.idx);
        conn[{v[0], v[2]}].push_back(f.idx);
    }

    DualGraph g = make_dual(face, conn);

    // remove connection between
    for (auto& kv: conn)
    {
        if (kv.second.size() > 3)
        {
            for (std::size_t i = 0; i < kv.second.size() - 1; ++i)
            {
                for (std::size_t j = i + 1; j < kv.second.size(); ++j)
                {
                    g.remove_edge(kv.second[i], kv.second[j]);
                }
            }
        }
    }

    // rm_iregular_node(&g);
    rm_debris(&g);
    std::vector<std::size_t> newf = g.allnode();
    std::set<std::size_t> newv{};
    for (auto& fi: newf)
    {
        newv.insert(face_org[fi].vertex.begin(), face_org[fi].vertex.end());
    }

    std::vector<std::size_t> vid_map(vertex_org.size());
    auto it = newv.begin();
    for (std::size_t i = 0; i < newv.size(); ++i)
    {
        vid_map[*it] = i;
        ++it;
    }

    std::vector<std::array<double, 3>> retv{};
    retv.reserve(newv.size());
    std::vector<std::array<std::size_t, 3>> retf{};
    retf.reserve(newf.size());
    for (const auto vi: newv)
    {
        const auto& v = vertex_org[vi];
        retv.push_back({v.x, v.y, v.z});
    }
    for (const auto& fi: newf)
    {
        const auto& v = face_org[fi].vertex;
        retf.push_back({vid_map[v[0]], vid_map[v[1]], vid_map[v[2]]});
    }

    return {retv, retf};
}

} // end of namespace post_acvd
