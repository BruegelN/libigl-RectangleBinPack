#include <array>
#include <vector>
#include <set>
#include <cmath> // std::floor, std::ceil

#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/facet_components.h>
#include <igl/remove_unreferenced.h>
#include <igl/edges.h>
#include <igl/cut_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/lscm.h>
#include <igl/combine.h>
#include <igl/avg_edge_length.h>

#include "Rect.h"
#include "GuillotineBinPack.h"



Eigen::MatrixXd V, V_uv;
Eigen::MatrixXi F, F_uv;
std::vector<std::vector<int>> cuts;
std::vector<Eigen::MatrixXd> vecV, vecV_uv;
std::vector<Eigen::MatrixXi> vecF, vecF_uv;

void initialize_segmentation( const Eigen::MatrixXd &V,const Eigen::MatrixXi &F,
                              std::vector<std::vector<int>> &cuts);

void cut(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<std::vector<int>> &cuts,
         std::vector<Eigen::MatrixXd> &vecV, std::vector<Eigen::MatrixXi> &newF);

void LSCM(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
          Eigen::MatrixXd &V_uv, Eigen::MatrixXi &F_uv);

void mergeUV(const std::vector<Eigen::MatrixXd> &vecV_uv,const std::vector<Eigen::MatrixXi> &vecF_uv,
            Eigen::MatrixXd &V_uv, Eigen::MatrixXi &F_uv);

void scale( std::vector<Eigen::MatrixXd> &vecV, std::vector<Eigen::MatrixXi> vecF,
            std::vector<Eigen::MatrixXd> &vecV_uv, std::vector<Eigen::MatrixXi> &vecF_uv);

void pack(std::vector<Eigen::MatrixXd> &vecV);

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;

    if (key == '0')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V_uv, F_uv);
        viewer.core().align_camera_center(V_uv, F_uv);
    }
    else if (key == ' ')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.core().align_camera_center(V, F);
    }
    else if((key-'1') >= 0 && (key-'1') < vecV_uv.size())
    {
        int i = key-'1'; // key '1' should be 1-th element(vec.at(0)!)
        viewer.data().clear();
        viewer.data().set_mesh(vecV_uv.at(i), vecF_uv.at(i));
        viewer.core().align_camera_center(vecV_uv.at(i), vecF_uv.at(i));
    }
    else
        return false;
    return true;
}


int main(int argc, char *argv[])
{
    if(argc <= 1)
    {
        std::cerr << "please provide the path to a triangle mesh as argument\n";
        std::cerr << argv[0] <<" ~/foo/data/mesh.obj\n";
        return EXIT_FAILURE;
    }

    igl::read_triangle_mesh(argv[1],V,F);

    initialize_segmentation(V,F,cuts);

    cut(V, F, cuts, vecV, vecF);

    vecV_uv.clear();
    vecF_uv.clear();
    for(size_t i = 0; i < vecV.size(); i++)
    {
        // parameterize per patch
        Eigen::MatrixXd v = vecV.at(i), v_uv;
        Eigen::MatrixXi f = vecF.at(i), f_uv;

        LSCM(v,f, v_uv, f_uv);

        vecV_uv.push_back(v_uv);
        vecF_uv.push_back(f_uv);
    }

    scale(vecV,vecF,vecV_uv,vecF_uv);

    // collect all pre- centers, pack, and then move em
    pack(vecV_uv);

    igl::combine(vecV_uv, vecF_uv, V_uv, F_uv);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_down = &key_down;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();
}

void initialize_segmentation(const Eigen::MatrixXd& V,const Eigen::MatrixXi& F,
                             std::vector<std::vector<int>> & cuts)
{
    Eigen::VectorXd X;
    Eigen::MatrixXi E;
    igl::edges(F,E);
    X = V.col(0);
    for(size_t i = 0; i < E.rows(); i++)
    {
        if( (X(E(i,0)) < X.mean() && X(E(i,1)) > X.mean()) ||
            (X(E(i,0)) > X.mean() && X(E(i,1)) < X.mean()))
        {
            std::vector<int> c;
            c.push_back(E(i,0));
            c.push_back(E(i,1));
            cuts.push_back(c);
        }
    }
    X = V.col(1);
    for(size_t i = 0; i < E.rows(); i++)
    {
        if( (X(E(i,0)) < X.mean() && X(E(i,1)) > X.mean()) ||
            (X(E(i,0)) > X.mean() && X(E(i,1)) < X.mean()))
        {
            std::vector<int> c;
            c.push_back(E(i,0));
            c.push_back(E(i,1));
            cuts.push_back(c);
        }
    }
    X = V.col(2);
    for(size_t i = 0; i < E.rows(); i++)
    {
        if( (X(E(i,0)) < X.mean() && X(E(i,1)) > X.mean()) ||
            (X(E(i,0)) > X.mean() && X(E(i,1)) < X.mean()))
        {
            std::vector<int> c;
            c.push_back(E(i,0));
            c.push_back(E(i,1));
            cuts.push_back(c);
        }
    }

}

void split(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
           std::vector<Eigen::MatrixXd> &Vs, std::vector<Eigen::MatrixXi> &Fs)
{
    Vs.clear();
    Fs.clear();

    Eigen::VectorXi C;
    igl::facet_components(F,C);
    std::vector<Eigen::MatrixXi> vec;
    for(int i = 0; i <= C.maxCoeff(); i++)
    {
        vec.emplace_back(1,3);
    }

    for(int rowIndex = 0; rowIndex < F.rows(); rowIndex++)
    {
        int comp = C(rowIndex);
        auto face = F.row(rowIndex);

        vec.at(comp).conservativeResize(vec.at(comp).rows()+1, vec.at(comp).cols());
        vec.at(comp).row(vec.at(comp).rows()-1) = face;
    }

    for(int i = 0; i <= C.maxCoeff(); i++)
    {
        Eigen::MatrixXi f = vec.at(i).bottomRows(vec.at(i).rows()-1);
        Eigen::MatrixXd v = V;
        Eigen::VectorXi I;
        igl::remove_unreferenced(Eigen::MatrixXd(v),Eigen::MatrixXi(f),v,f,I);
        Vs.push_back(v);
        Fs.push_back(f);
    }
    assert(Vs.size() == Fs.size());
}

void cut(const Eigen::MatrixXd& V, const Eigen::MatrixXi &F, const std::vector<std::vector<int>> &cuts,
         std::vector<Eigen::MatrixXd> &Vs, std::vector<Eigen::MatrixXi> &Fs)
{
    const size_t num_faces = F.rows();
    Eigen::MatrixXi cut_mask(num_faces, 3);
    std::set<std::array<int, 2>> cut_edges;
    for (const auto& cut : cuts)
    {
        const size_t cut_len = cut.size();
        for (size_t i=0; i<cut_len-1; i++)
        {
            std::array<int, 2> e{cut[i], cut[i+1]};
            if (e[0] > e[1])
            {
                std::swap(e[0], e[1]);
            }
            cut_edges.insert(e);
        }
    }

    cut_mask.setZero();
    for (size_t i=0; i<num_faces; i++)
    {
        std::array<int, 2> e0{F(i, 0), F(i, 1)};
        std::array<int, 2> e1{F(i, 1), F(i, 2)};
        std::array<int, 2> e2{F(i, 2), F(i, 0)};
        if (e0[0] > e0[1]) std::swap(e0[0], e0[1]);
        if (e1[0] > e1[1]) std::swap(e1[0], e1[1]);
        if (e2[0] > e2[1]) std::swap(e2[0], e2[1]);

        if (cut_edges.find(e0) != cut_edges.end()) {
            cut_mask(i, 0) = 1;
        }
        if (cut_edges.find(e1) != cut_edges.end()) {
            cut_mask(i, 1) = 1;
        }
        if (cut_edges.find(e2) != cut_edges.end()) {
            cut_mask(i, 2) = 1;
        }
    }

    Eigen::MatrixXd V_new;
    Eigen::MatrixXi F_new;
    igl::cut_mesh(V,F,cut_mask, V_new, F_new);
    split(V_new, F_new, Vs, Fs);

}


void LSCM(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
          Eigen::MatrixXd &V_uv, Eigen::MatrixXi &F_uv)
{
        std::vector<std::vector<int>> boundary_loop_list;
        igl::boundary_loop(F, boundary_loop_list);

        int loop_length = 0;
        std::vector<int> largest_boundary_loop;
        for (auto &boundary_loop : boundary_loop_list)
        {
            // assume that more vertices along a boundary correspond to longer boundary
            if (loop_length < boundary_loop.size())
            {
                loop_length = boundary_loop.size();
                largest_boundary_loop = boundary_loop;
            }
        }

        Eigen::MatrixXd bnd_uv;
        Eigen::Map<Eigen::VectorXi> bnd(largest_boundary_loop.data(), largest_boundary_loop.size());

        // Fix two points on the boundary
        Eigen::VectorXi b(2, 1);
        b(0) = bnd(0);
        b(1) = bnd(round(bnd.size() / 2));
        double distance = Eigen::VectorXd((V.row(b(1))).array()-(V.row(b(0))).array()).norm();
        Eigen::MatrixXd bc(2,2);
        bc<< 0,0,0,distance;


        // LSCM parametrization
        igl::lscm(V, F, b, bc, V_uv);

        F_uv = F;
}

void scale( std::vector<Eigen::MatrixXd> &vecV, std::vector<Eigen::MatrixXi> vecF,
            std::vector<Eigen::MatrixXd> &vecV_uv, std::vector<Eigen::MatrixXi> &vecF_uv)
{
    if(vecV.size() != vecV_uv.size())
    {
        std::cerr << "vecV_uv.size() != vecF.size()" << std::endl;
    }

    for(size_t i = 0; i < vecV.size(); i++)
    {
        double avg = igl::avg_edge_length(vecV.at(i), vecF.at(i));
        double avg_uv = igl::avg_edge_length(vecV_uv.at(i), vecF_uv.at(i));
        vecV_uv.at(i) = (vecV_uv.at(i).array()*avg/avg_uv);
    }
}

void pack(std::vector<Eigen::MatrixXd> &vecV_uv)
{
    const size_t nPatches = vecV_uv.size();
    constexpr double scale_to_int = 1000.0;
    if (100 < nPatches)
    {
        std::cerr << "Try to pack more than 100 patches, it might be slow!" << std::endl;
    }

    std::vector<Eigen::RowVector2d> centersPre(nPatches);
    std::vector<Eigen::RowVector2d> centersPost(nPatches);
    std::vector<Eigen::RowVector2d> minCoeffs(nPatches);
    std::vector<Eigen::RowVector2d> maxCoeffs(nPatches);
    std::vector<rbp::Rect> rects(nPatches);
    std::vector<bool> vecFlipped(nPatches);

    std::vector<Eigen::MatrixXi> binPackVs(nPatches);
    double area = 0.0;
    for(size_t i = 0;  i < nPatches; i++)
    {
        Eigen::MatrixXd v = vecV_uv.at(i);
        Eigen::RowVectorXd min = v.colwise().minCoeff();
        // move min to (0,0) such that all values are >= 0
        v = v.rowwise() - min;
        // largest coordinates can now be used at edge length for BBox
        Eigen::RowVectorXd max = v.colwise().maxCoeff();
        // + 5% bounardy
        const double width = max(0) * 1.05;
        const double height = max(1)* 1.05;
        // accumulate area for a first guess for the total area required
        area += width*height;
        // make  BBox centered at (0,0)
        Eigen::RowVector2d center;
        center << (width/2.0), (height/2.0);
        centersPre.at(i) = center;
        vecV_uv.at(i) = v.rowwise() - center;
        // convert to int
        int iWidth = static_cast<int>(std::ceil(width*scale_to_int));
        int iHeight = static_cast<int>(std::ceil(height*scale_to_int));
        rects.at(i) = rbp::Rect{0, 0, iWidth, iHeight};
    }

    {
        int box_side = static_cast<int>(std::ceil(std::sqrt(area * scale_to_int * 1.02)));
        bool failed;
        size_t num_faild = 0;
        constexpr size_t max_tries = 100;
        do
        {
            // reset everytime, in case failed befores
            failed  = false;
            //vecFlipped.clear();

            rbp::GuillotineBinPack packer(box_side, box_side);
            for (size_t i = 0; i < nPatches; i++) {
                // try to insert one rectangle after the other
                rbp::Rect curr_rect = rects.at(i);
                rbp::Rect result = packer.Insert(curr_rect.width, curr_rect.height, true,
                                                 rbp::GuillotineBinPack::FreeRectChoiceHeuristic::RectBestAreaFit,
                                                 rbp::GuillotineBinPack::GuillotineSplitHeuristic::SplitMinimizeArea);
                if (0 == result.width)
                {
                    ++num_faild;
                    failed = true;
                    break;
                }

                // need to rotate later on
                vecFlipped.at(i) = (curr_rect.height == result.width);
                // from int back to original scale and safe for later
                Eigen::RowVector2d center;
                center << (result.x+result.width/2), (result.y+result.height/2);
                center = center.array()/scale_to_int;
                centersPost.at(i) = center;
            }
            // increase in 10% steps, if it failed
            box_side *= 1.1;
        } while (failed && num_faild < max_tries);
        if(num_faild >= max_tries)
        {
            std::cerr << "Reached max number of tries to pack UV!" << std::endl;
        }
    }

    for(size_t i = 0; i < nPatches; i++)
    {
        Eigen::MatrixXd v = vecV_uv.at(i);

        if(true == vecFlipped.at(i))
        {
            const Eigen::MatrixXd rotation =(Eigen::MatrixXd(2,2)<<
                0.0, 1.0,
                -1.0, 0.0).finished();
            v = (v*rotation).eval();  // rotation.T not required bcs BBox is centered at (0,0)
        }
        // just move the element at the right place
        v = (v.rowwise() + centersPost.at(i)).eval();
        vecV_uv.at(i) = v;
    }
}

