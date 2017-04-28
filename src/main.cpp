#include "Mesh.hpp"
#include "Viewer.hpp"
#undef Success  
#include "Eigen/Sparse"
#include <cstdlib>
#include "Eigen/SparseLU"
#include <vector>
#include "DGP/Plane3.hpp"

// #include "Eigen/UmfPackSupport"
// #include "Eigen/umfpack.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;

void calcuateS(Mesh & source_mesh, Mesh & deformed_source_mesh, vector<Matrix3> & S);
void calcuateA_c(vector<pair<long,MeshFace *> > &M,vector<Matrix3> &S, Mesh & target_mesh, Eigen::SparseMatrix<double,Eigen::RowMajor> & A, Eigen::SparseVector<double> & c);
void calculateM(string correspondence_path, Mesh & source_mesh, Mesh & target_mesh, vector<pair<long,MeshFace *> > &M);

class MatrixReplacement;
using Eigen::SparseMatrix;
namespace Eigen {
namespace internal {
  // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<MatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};
}
}
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };
  Index rows() const { return mp_mat->rows(); }
  Index cols() const { return mp_mat->cols(); }
  template<typename Rhs>
  Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
  // Custom API:
  MatrixReplacement() : mp_mat(0) {}
  void attachMyMatrix(const SparseMatrix<double> &mat) {
    mp_mat = &mat;
  }
  const SparseMatrix<double> my_matrix() const { return *mp_mat; }
private:
  const SparseMatrix<double> *mp_mat;
};
// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
  template<typename Rhs>
  struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixReplacement,Rhs,generic_product_impl<MatrixReplacement,Rhs> >
  {
    typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      assert(alpha==Scalar(1) && "scaling is not implemented");
      // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
      // but let's do something fancier (and less efficient):
      for(Index i=0; i<lhs.cols(); ++i)
        dst += rhs(i) * lhs.my_matrix().col(i);
    }
  };
}
}


int
usage(int argc, char * argv[])
{
  DGP_CONSOLE << "";
  DGP_CONSOLE << "Usage: " << argv[0] << " <source-in> <source-deformed> <target-in> <source_target_correspondence>";
  DGP_CONSOLE << "";

  return -1;
}

int
main(int argc, char * argv[])
{
  if (argc < 2)
    return usage(argc, argv);

  if (argc == 3){
    std::string path = argv[1];
    std::string out_path = argv[2];
 
   Mesh mesh ;
    if (!mesh.load(path)){
      return -1;
    }
    std::ofstream out(out_path.c_str(), std::ios::binary);
    if (!out)
    {
      DGP_ERROR << "Could not open '" << out_path << "' for writing";
      return -1;
    }

   for(auto it = mesh.facesBegin(); it!=mesh.facesEnd(); ++it){
      auto vit = (*it).verticesBegin(); 
      Vector3 v0 = (*vit)->getPosition();
      vit++;
      Vector3 v1 = (*vit)->getPosition();
      vit++;
      Vector3 v2 = (*vit)->getPosition(); 
      Vector3 point  = (v0+v1+v2)/3;
      Vector3 normal = (*it).getNormal();
      out << point[0] << ' ' << point[1] << ' ' << point[2] << ' ' <<normal[0] << ' ' << normal[1] << ' ' << normal[2] << '\n';
    }
  return 0;
  }

  std::string source_path = argv[1];
  std::string source_deformed_path = argv[2];
  std::string target_path = argv[3];
  std::string correspondence_path = argv[4];

  Mesh source_mesh, deformed_source_mesh, target_mesh, deformed_target_mesh;
  
  if (!source_mesh.load(source_path) || !deformed_source_mesh.load(source_deformed_path) || !target_mesh.load(target_path) || !deformed_target_mesh.load(target_path))
    return -1;

  DGP_CONSOLE << "Read mesh '" << source_mesh.getName() << "' with " << source_mesh.numVertices() << " vertices, " << source_mesh.numEdges()
              << " edges and " << source_mesh.numFaces() << " faces from " << source_path;
  DGP_CONSOLE << "Read mesh '" << deformed_source_mesh.getName() << "' with " << deformed_source_mesh.numVertices() << " vertices, " << deformed_source_mesh.numEdges()
              << " edges and " << deformed_source_mesh.numFaces() << " faces from " << source_deformed_path;
  DGP_CONSOLE << "Read mesh '" << target_mesh.getName() << "' with " << target_mesh.numVertices() << " vertices, " << target_mesh.numEdges()
              << " edges and " << target_mesh.numFaces() << " faces from " << target_path;


  //setting id of meshvertex
  long li = 0;
  for(auto it = source_mesh.verticesBegin(); it!=source_mesh.verticesEnd(); ++it, li++)it->id = li;

  li = 0;
  for(auto it = deformed_source_mesh.verticesBegin(); it!=deformed_source_mesh.verticesEnd(); ++it, li++)it->id = li;

  li = 0;
  for(auto it = target_mesh.verticesBegin(); it!=target_mesh.verticesEnd(); ++it, li++)it->id = li;

  li = 0;
  for(auto it = deformed_target_mesh.verticesBegin(); it!=deformed_target_mesh.verticesEnd(); ++it, li++)it->id = li;


  //for(auto it = source_mesh.verticesBegin(); it!=source_mesh.verticesEnd(); ++it)DGP_CONSOLE << it->id;

  //for(auto it = deformed_source_mesh.verticesBegin(); it!=deformed_source_mesh.verticesEnd(); ++it)DGP_CONSOLE << it->id;

  //for(auto it = target_mesh.verticesBegin(); it!=target_mesh.verticesEnd(); ++it)DGP_CONSOLE << it->id;

  //for(auto it = deformed_target_mesh.verticesBegin(); it!=deformed_target_mesh.verticesEnd(); ++it)DGP_CONSOLE << it->id;

  // if (target_num_faces >= 0 && mesh.numFaces() > target_num_faces)
  // {
  //   mesh.updateBounds();
  // }

  // if (!out_path.empty())
  // {
  //   if (!mesh.save(out_path))
  //     return -1;

  //   DGP_CONSOLE << "Saved mesh to " << out_path;
  // }

  // source_mesh.check(deformed_source_mesh);			//check source and deformed dimensions and also check for triangular mesh

  vector<Matrix3> S;								//source deformations
  vector<pair<long,MeshFace *> > M;
  // M.resize(10);

  calculateM(correspondence_path,source_mesh,target_mesh,M);

  // auto it = target_mesh.facesBegin();
  
  // for(int i=0;it!=target_mesh.facesEnd();i++,it++)
  // {
  // 	M.push_back(pair<long,MeshFace*>(i,&(*it)));
  // }

  //DGP_CONSOLE<<"1\n";
  Eigen::SparseMatrix<double,Eigen::RowMajor> A(9*M.size(),3*(target_mesh.numVertices()+target_mesh.numFaces())); 
  //DGP_CONSOLE<<"2\n";
  Eigen::SparseVector<double> c(9*M.size());
  //DGP_CONSOLE<<"3\n";
  calcuateS(source_mesh,deformed_source_mesh,S);	//calcuate S
  //DGP_CONSOLE<<"4\n";
  calcuateA_c(M,S,target_mesh,A,c);
  //DGP_CONSOLE<<"5\n";

  //DGP_CONSOLE << "A rows\t" << A.rows() << "\n";
  //DGP_CONSOLE << "A cols\t" << A.cols() << "\n";
  //DGP_CONSOLE << "c rows\t" << c.rows() << "\n";
  //DGP_CONSOLE << "c cols\t" << c.cols() << "\n";

  //calculating A(transpose)*A
  // Eigen::SparseMatrix<double,Eigen::RowMajor> AtA(3*(target_mesh.numVertices()+target_mesh.numFaces()), 3*(target_mesh.numVertices()+target_mesh.numFaces())); 
  // AtA = A.transpose() * A;

  //DGP_CONSOLE << "AtA rows\t" << AtA.rows() << "\n";
  //DGP_CONSOLE << "AtA cols\t" << AtA.cols() << "\n";

  //calculating Atc
  // Eigen::SparseVector<double> Atc(3*(target_mesh.numVertices()+target_mesh.numFaces()));
  // Atc = A.transpose() * c;

  //DGP_CONSOLE << "Atc rows\t" << Atc.rows() << "\n";
  //DGP_CONSOLE << "Atc cols\t" << Atc.cols() << "\n";

  // Eigen::SparseLU<SparseMatrix<double,Eigen::RowMajor>, Eigen::COLAMDOrdering<Eigen::Index> > solver;

  // Eigen::SparseMatrix<double,Eigen::RowMajor> L(3*(target_mesh.numVertices()+target_mesh.numFaces()), 3*(target_mesh.numVertices()+target_mesh.numFaces())); 
  // Eigen::SparseMatrix<double,Eigen::RowMajor> U(3*(target_mesh.numVertices()+target_mesh.numFaces()), 3*(target_mesh.numVertices()+target_mesh.numFaces())); 

  // L = AtA.triangularView<Eigen::StrictlyLower>();
  // U = AtA.triangularView<Eigen::Upper>();

  //DGP_CONSOLE << "L rows\t" << L.rows() << "\n";
  //DGP_CONSOLE << "L cols\t" << L.cols() << "\n";

  //DGP_CONSOLE << "U rows\t" << U.rows() << "\n";
  //DGP_CONSOLE << "U cols\t" << U.cols() << "\n";

  //DGP_CONSOLE << "A != 0\t" << A.nonZeros() << "\n";
  //DGP_CONSOLE << "c != 0\t" << c.nonZeros() << "\n";

  //DGP_CONSOLE << "AtA != 0\t" << AtA.nonZeros() << "\n";
  //DGP_CONSOLE << "Atc != 0\t" << Atc.nonZeros() << "\n";

  //Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::UmfPack> lu_of_A(AtA);

  //VectorXd x(AtA.rows());
  //Eigen::matrixL().solve(Atc);

  //x = Eigen::SparseLU<SparseMatrix<double, Eigen::RowMajor>, Eigen::COLAMDOrdering<int> >(AtA).solve(Atc);

  //AtA.finalize();
  //Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::UmfPack> lu_of_A(AtA);

  //solver.analyzePattern(AtA);
  //solver.factorize(AtA);
  //Eigen::VectorXd x;
  //x = solver.solve(Atc);

  /*MatrixReplacement AA;
  AA.attachMyMatrix(A);

  Eigen::VectorXd x;

  Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
  cg.compute(AA);
  x = cg.solve(c);
  DGP_CONSOLE << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << "\n";
*/
  //DGP_CONSOLE<<A;



  Viewer viewer1, viewer2, viewer3;
  viewer1.setObject(&source_mesh);
  viewer1.launch(argc, argv);

  return 0;
}


void calcuateS(Mesh & source_mesh, Mesh & deformed_source_mesh, vector<Matrix3> & S)
{
	for(Mesh::FaceIterator fit = source_mesh.facesBegin(), dfit = deformed_source_mesh.facesBegin(); fit!= source_mesh.facesEnd() && dfit!=deformed_source_mesh.facesEnd(); fit++, dfit++){
		list<Mesh::Vertex *>::iterator  vit = fit->verticesBegin();
		Vector3 v1 = (*vit)->getPosition();
		vit++;
		Vector3 v2 = (*vit)->getPosition();
		vit++;
		Vector3 v3 = (*vit)->getPosition();
		Vector3 v4 = v1 + ((v2 - v1).cross((v3 - v1)))/((v2 - v1).cross(v3 - v1)).length();

		list<Mesh::Vertex *>::iterator dvit = dfit->verticesBegin();
		Vector3 dv1 = (*dvit)->getPosition();
		dvit++;
		Vector3 dv2 = (*dvit)->getPosition();
		dvit++;
		Vector3 dv3 = (*dvit)->getPosition();
		Vector3 dv4 = dv1 + ((dv2 - dv1).cross((dv3 - dv1)))/((dv2 - dv1).cross(dv3 - dv1)).length();

		Matrix3  V = Matrix3::fromColumns(v2-v1,v3-v1,v4-v1);  
		Matrix3  dV = Matrix3::fromColumns(dv2-dv1,dv3-dv1,dv4-dv1);  

		if(V.determinant()<1e-10){
			DGP_CONSOLE<<"determinant error"<<std::endl;
			exit(0);
		}
		S.push_back(dV*V.inverse());
	}
}

void calcuateA_c(vector<pair<long,MeshFace*> > &M,vector<Matrix3> &S, Mesh & target_mesh, Eigen::SparseMatrix<double,Eigen::RowMajor> & A, Eigen::SparseVector<double> & c)
{

  ofstream a_file;
  ofstream c_file;

  a_file.open("A.txt");
  c_file.open("c.txt");

  a_file << A.rows() << " " << A.cols() << "\n";
  c_file << c.rows() << " " << c.cols() << "\n"; 

	for(unsigned long i=0; i < M.size(); i++)
	{
		list<Mesh::Vertex *>::iterator  vit = M[i].second->verticesBegin();
		Vector3 v1 = (*vit)->getPosition();
		MeshVertex *vert1 = (*vit);
		vit++;
		Vector3 v2 = (*vit)->getPosition();
		MeshVertex *vert2 = (*vit);
		vit++;
		Vector3 v3 = (*vit)->getPosition();
		MeshVertex *vert3 = (*vit);
		Vector3 v4 = v1 + ((v2 - v1).cross((v3 - v1)))/((v2 - v1).cross(v3 - v1)).length();
		Matrix3 V = Matrix3::fromColumns(v2-v1,v3-v1,v4-v1);  

		if(V.determinant()<1e-10)
		{
			DGP_CONSOLE<<"determinant error"<<std::endl;
			exit(0);
		}

		//DGP_CONSOLE<<"AC1\n";
    // v inverse for target
		V.invert();
		//DGP_CONSOLE<<"AC2\n";

    // DGP_CONSOLE << vert1->id << " " << vert2->id << " " << vert3->id << "\n";

		A.reserve(Eigen::VectorXi::Constant(A.rows(),12));
		for(int j=0;j<3;j++)
		{
			for(int k=0;k<3;k++)
			{
				// cout<<i<<" "<<j<<" "<<k<<endl;
				// [x1 x2 x3] * for all vertices then [x y z] of v4 corresponding to all triangles
				
        // setting coefficients of c
        c.coeffRef(i*9+j*3+k) = S[M[i].first](j,k);

        c_file << i*9+j*3+k << " " << S[M[i].first](j,k) << "\n";

        // setting coefficients of A
				A.insert(i*9+j*3+k,vert1->id*3+j) = -(V(0,k)+V(1,k)+V(2,k));
        a_file << i*9+j*3+k << " " << vert1->id*3+j << " " << -(V(0,k)+V(1,k)+V(2,k)) << "\n";

				A.insert(i*9+j*3+k,vert2->id*3+j) = V(0,k);
        a_file <<  i*9+j*3+k << " " << vert2->id*3+j << " " << V(0,k) << "\n";

				A.insert(i*9+j*3+k,vert3->id*3+j) = V(1,k);
        a_file << i*9+j*3+k << " " << vert3->id*3+j << " " << V(1,k) << "\n";

				A.insert(i*9+j*3+k,target_mesh.numVertices()*3 + M[i].second->id*3 + j) = V(2,k);
        a_file << i*9+j*3+k << " " << target_mesh.numVertices()*3 + M[i].second->id*3 + j << " " << V(2,k) << "\n";

			}
		}

	}

  a_file.close();
  c_file.close();

}


void calculateM(string correspondence_path, Mesh & source_mesh, Mesh & target_mesh, vector<pair<long,MeshFace *> > &M)
{
    std::ifstream in(correspondence_path.c_str());
    if (!in)
    {
      DGP_ERROR << "Could not open correspondence map" << correspondence_path << "' for reading";
      return;
    } 

  long np;
    if (!(in >>np))
    {
      DGP_ERROR << "Could not read number of correspondence points from file '" << correspondence_path << '\'';
      return;
    }
    // clear();

    vector<Vector3> source_vertices;
    vector<Vector3> target_vertices;
    Vector3 s;
    Vector3 t;
    for (long i = 0; i < np; ++i)
    {
      if (!(in >> t[0] >> t[1] >> t[2] >> s[0] >> s[1] >> s[2]))
      {
        DGP_ERROR << "Could not read correspondence points ";
      }
      source_vertices.push_back(s);
      target_vertices.push_back(t);
    }
    // cout<<"source"<<endl;;
    // for(long i =0 ;i < source_vertices.size();i++){
    //   cout<<source_vertices[i][0]<<" "<<source_vertices[i][1]<<" "<<source_vertices[i][2]<<endl;
    // }
    // cout<<"target"<<endl;
    // for(long i =0 ;i < target_vertices.size();i++){
    //   cout<<target_vertices[i][0]<<" "<<target_vertices[i][1]<<" "<<target_vertices[i][2]<<endl;
    // }
    if(target_mesh.numFaces() != np){
        DGP_ERROR<<"Number of correspondence points is not equal to number of faces of target mesh";
        return ;
    }
    auto tit = target_mesh.facesBegin();
    for(long i=0;i<np;i++,tit++){
      long source_index = 0;
      MeshFace * target_face = &(*tit) ;

      double min_dist = numeric_limits<double>::max();

      for(auto it = source_mesh.facesBegin(); it!= source_mesh.facesEnd(); it++){
          auto vit = (*it).verticesBegin(); 
          Vector3 v0 = (*vit)->getPosition();
          vit++;
          Vector3 v1 = (*vit)->getPosition();
          vit++;
          Vector3 v2 = (*vit)->getPosition();
          double temp_dist = (v0 - source_vertices[i]).length() + (v1 - source_vertices[i]).length() + (v2 - source_vertices[i]).length();
          if(temp_dist < min_dist){
            min_dist = temp_dist;
            source_index = (*it).id;
          }
        }     

      M.push_back(pair<long,MeshFace*>(source_index,target_face));
    }   
}