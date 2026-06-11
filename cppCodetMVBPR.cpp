//c++ Script tratment selection via MultiView Bayesian Profile Regression
//author: Lorenzo Moni, Silvia Liverani, Alberto Cassese, Francesco Claudio Stingo  



#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp.h>
#include <chrono>  // Include chrono for timing
#include <stack>
#include <cmath>
#include <optional>


using namespace Rcpp;
using namespace std;
//[[Rcpp::depends(RcppArmadillo, RcppDist)]]

//simple hash function for tuple. This function can be optimize to have a perfect hash  
struct tuple_hash {
  arma::uword operator()(const std::tuple<arma::uword, arma::uword, arma::uword>& t) const {
    arma::uword h1 = std::get<0>(t);
    arma::uword h2 = std::get<1>(t);
    arma::uword h3 = std::get<2>(t);
    
    // Combine the hashes using a robust method
    arma::uword seed = 0;
    seed ^= h1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    
    return seed;
  }
};

struct tuple_hashperfect {
  // Bit widths for each component
  static constexpr arma::uword bv = 10; // v: 0..1023
  static constexpr arma::uword bk = 27; // k: 0..134,217,727
  static constexpr arma::uword bd = 27; // d: 0..134,217,727
  
  arma::uword operator()(const std::tuple<arma::uword, arma::uword, arma::uword>& t) const {
    arma::uword v = std::get<0>(t);
    arma::uword k = std::get<1>(t);
    arma::uword d = std::get<2>(t);
    
    // Bit-packing into 64-bit integer
    // layout: [ d(27) | k(27) | v(10) ]
    return (d << (bk + bv)) | (k << bv) | v;
  }
};


//-----------Print function

void printtuplelist(std::unordered_map<std::tuple<arma::uword, arma::uword, arma::uword>,
                    arma::Col<double>,
                    tuple_hash> SS){
  cout << "+++START PRINT+++"<<endl;
  for (auto s : SS) {
    // Unpack the tuple and print its elements
    auto key = s.first;
    cout << "S_vkd: ("
    << std::get<0>(key) << ", "
    << std::get<1>(key) << ", "
    << std::get<2>(key) << ") "
    << "S_vk: " << s.second.t() << endl;  // Transpose for better formatting
  }
  cout << "+++END PRINT+++"<<endl;
}


template<typename T>
  void printUnordlist(const std::unordered_map<arma::uword, T>& um){
    
    for (auto s : um){
      cout<< "key: " << s.first  << " value: "<<s.second <<endl;
      
    }
  }
void printUnordMapVec(const std::unordered_map<arma::uword, arma::Col<double>>& um) {
  for (const auto& [key, vec] : um) {
    std::cout << "key: " << key << " value: (";
    
    // Print the vector elements within parentheses, separated by commas
    for (arma::uword i = 0; i < vec.n_elem; ++i) {
      std::cout << vec(i);
      if (i != vec.n_elem - 1) {
        std::cout << ", ";  // Add a comma between values
      }
    }
    
    std::cout << ")\n";  // Close the parentheses and move to the next line
  }
}


void printSet(const std::unordered_set<arma::uword>& s) {
  std::cout << "{ ";
  for (const auto& elem : s) {
    std::cout << elem << " ";
  }
  std::cout << "}" << std::endl;
}



void printContainerNewIndices(std::vector<std::stack<arma::uword>> c) {
  std::cout  <<endl;
  
  std::cout<< "print ContainerNewIndices "  <<endl;
  
  auto Nvector=c.size();
  
  for(arma::uword v=0; v<Nvector;v++){
    std::cout << "View: "<< v+1 << ": ";
    
    auto temp_stack1 = c.at(v);
    while (!temp_stack1.empty()) {
      std::cout << temp_stack1.top() << " ";
      temp_stack1.pop();
    }
    std::cout  <<endl;
    
  }
  std::cout  <<endl;
  
}

//--------
  class Ymodel{
    public:
      //hyperparameter tstudent theta_k,r r=1,...,R_y-1
      arma::Col<double> theta0, scale, df;//, NorConst1;
      arma::Col<double> inv_dfscale2,  NegdfPlusOnedivTwo;
      double logNorConst;
      arma::uword Ry, RylessTwo, RylessOne;
      //note that inv_dfscale, NegdfPlusOnedivTwo are vector containing the values
      //for each(independent) theta_k,r .. summing at the end will give us the complete log prior
      
      
      Ymodel(Rcpp::List list, arma::uword ryless1)
      {
        
        if (list.containsElementNamed("HParthetak")) {
          // should be list of list
          Rcpp::List listHYtheta=list["HParthetak"];
          //theta0
          if(listHYtheta.containsElementNamed("mean")){
            theta0=as<arma::Col<double>>(listHYtheta["mean"]);
          }else{
            throw std::runtime_error("hyperparameter mean for theta_k's not given");
          }
          //scale
          if(listHYtheta.containsElementNamed("scale")){
            scale=as<arma::Col<double>>(listHYtheta["scale"]);
          }else{
            throw std::runtime_error("hyperparameter scale for theta_k's not given");
          }
          //df
          if(listHYtheta.containsElementNamed("df")){
            df=as<arma::Col<double>>(listHYtheta["df"]);
          }else{
            throw std::runtime_error("hyperparameter df for theta_k's not given");
          }
          
          //checks
          if(theta0.n_elem!=ryless1 || scale.n_elem!=ryless1||  df.n_elem!=ryless1){
            throw std::runtime_error("hyperparameter for theta_k: theta0 |or| df |or| scale.n_elem!=ryless1");
          }
          
          
          if(any(df < 0) || arma::any(df < 0)){
            throw std::runtime_error("hyperparameter for theta_k: df |or| scale <0");
          }
        }else{
          throw std::runtime_error("hyperparameters of theta_k's not given");
        }
        
        /*cout<< "theta " << theta0 <<endl;
        cout<< "scale " << scale <<endl;
        cout<< "df " << df <<endl;*/
          
          Ry=ryless1+1;
        RylessOne=ryless1;
        RylessTwo=ryless1-1;
        
        //pre-compute normalizing constant and inv_dfscale
        inv_dfscale2=1.0 / (df % arma::square(scale));
        NegdfPlusOnedivTwo= -0.5*(df+1);
        //NorConst= arma::accu(arma::lgamma(0.5*(df+1))-arma::lgamma(0.5*df)-0.5*inv_dfscale2);
        logNorConst=arma::accu( arma::lgamma(0.5*(df+1))-arma::lgamma(0.5*df)-
                                  0.5*arma::log(arma::datum::pi*df) - arma::log(scale));

        
        
        cout << "-> Initialization  yModel class........................................ OK" <<endl;
        
      }
      
      
      
      
      //return prior t-student computed in theta_kr
      // as prior we assume independent theta_kr~lst(theta_0r,s0r,df)
      double logPriorindtstud(const arma::Col<double>  & thetak){
        //LA COSTANTE NON MI DOVREBBE MAI SERVIRE MA PER ORA LA LASCIO
        
        
        arma::Col<double> VectLogKern=NegdfPlusOnedivTwo % arma::log(1+ (arma::square(thetak-theta0) % inv_dfscale2));

return (arma::accu(VectLogKern)+logNorConst);
      }
      
      //log likelihood of y
      //this method compute the log likelihood of y, that is
      // eta_i,y_i \propto exp(theta_z1,i + w_i'beta_r)
  double logLikcategoricalY(const arma::uword & yi,
                            const arma::Col<double> & thetak,
                            const arma::Col<double> & WBmat_t_allRi){
    //Vector of Rylessone+1 entries
    arma::Col<double>  ThetaplusWBallri(Ry);
    ThetaplusWBallri.zeros();

    //RylessTwo is the index to insert the fist Rylessone elements of thetak+WBmat_t_allRi
    // (W%*%Beta)^T [mat of R_y-i x n] of linear predictor of W'B_r by cols;
                               //WBmat_t_allRi = (W%*%Beta)^T .col(i) [col vect  W_i'B_1, W_i'B_2,...W_i'B_Ry-1 ]
    ThetaplusWBallri.subvec(0,RylessTwo)=thetak+WBmat_t_allRi;
    //not needed
    //ThetaplusWBallri(RylessOne)=0;

    //shift  ThetaplusWBallri-max
    ThetaplusWBallri-=ThetaplusWBallri.max();

    //sum_r exp(eta_r-max)
    double sumexp=arma::accu(arma::exp(ThetaplusWBallri));

    return ThetaplusWBallri(yi)-std::log(sumexp);
  }


  //sample from the the prior of theta_k, ie  G0_thetak a t-student
  arma::Mat<double> SampleG0thetak(arma::uword Msize=1){
    //it sample RY-1 components from the lst(theta0_r,s_r,df_r)
    //NOTE THE COMPONENTS ARE INDEPENDENTS
    //Tstud_lts(mu,scale,df)=mu+scale*[tstud(df)]=mu+scale*rnorm()/sqrt(rchisq(df)/df)
    //rchisq(df)/c =gamma(shape df/2, scale=2/c)

    //Matrix of iid from N(0,1)
    arma::Mat<double> Normstddraw=arma::randn(RylessOne, Msize);
    //Matrix of independent rchisq (the element r,m~rchisq(df=df_r))
    arma::Mat<double> InvSqrtSampChisqDfDivDf(RylessOne, Msize);
    for(arma::uword  ro=0; ro<InvSqrtSampChisqDfDivDf.n_rows;ro++){
      //sample gamma(df_r/2,2)=chisq(df)
      auto SampChisqDfDivDf= arma::randg<arma::Row<double>>(Msize, arma::distr_param(df(ro)/2,2.0));

      //compute sqrt(df/chisq) =1/sqrt SampChisqDfDivDf
      SampChisqDfDivDf=arma::sqrt(df(ro)/SampChisqDfDivDf);

      //store in InvSqrtSampChisqDfDivDf  SampChisqDfDivDf*scale
      InvSqrtSampChisqDfDivDf.row(ro)=SampChisqDfDivDf*scale(ro);
    }

    //element wise each col multiplication Norm*scale_r(sqrt(df_r/chisq_r))
    Normstddraw%=InvSqrtSampChisqDfDivDf;

    //add mu
    Normstddraw.each_col()+=theta0;
    return Normstddraw;

  }



  //Sample from the Y model: Y~cat(eta) with eta_r \propto exp(theta_kr+Wbeta_r)
  //Note that this modify the original matrix BUT is efficient since perform an inplace addition
  arma::Col<arma::uword> SampleResponse(arma::Mat<double>   & UnnormLogProbs){

    // Generate Gumbel noise: G = -log(-log(U)), U ~ Uniform(0,1) and
    //add to the matrix
    UnnormLogProbs+= -arma::log(-arma::log(arma::randu<arma::mat>(arma::size(UnnormLogProbs))));

    //note that this modify the  ORIGINAL UnnormLogProbs matrix (since UnnormLogProbs is a reference)  

    return arma::index_max(UnnormLogProbs, 0).t();

  }



};


//Class
//T student proposal
//propose a vector of scale.ncol
class MHProposal{
private:
  const double df;
  const double dim;

  arma::Mat<double> Scale, invScale;
  double logconst; //not useful
  arma::Col<double> Mzeros;



  //compute kernel of MV student t
  double logKernel(const arma::Col<double> &x,const arma::Col<double> &mu){
    //malanobis distance efficient

    arma::Col<double>  diff=x-mu;

    double MDdifdfpluss1= 1+arma::dot(diff, invScale*diff)/df;

    return -0.5*(dim+df)*std::log(MDdifdfpluss1);
  }


  //helper method init
  void SetScale(const arma::Mat<double> &scale_){
    Scale=scale_;
    invScale=arma::inv(Scale);

    logconst=std::lgamma((df+dim)*0.5)-std::lgamma(df*0.5)-(0.5*dim*std::log(arma::datum::pi*df));
    logconst+=-0.5*std::log(arma::det(Scale));


  }



public:
  //multivariate t student proposal constructor
  MHProposal(const arma::Mat<double> & scale_, double df_=5):
  df(df_), dim(static_cast<double>(scale_.n_cols)){
    //  invScale=arma::inv(Scale);

    // logconst=std::lgamma((df+dim)*0.5)-std::lgamma(df*0.5)-(0.5*dim*std::log(arma::datum::pi*df));
    //logconst+=-0.5*std::log(arma::det(Scale));

    SetScale(scale_);
    Mzeros.set_size(Scale.n_cols);
    Mzeros.zeros();
  }


  //change scale matrix
  void ChangeScaleProposalUsingCovMat(const arma::Mat<double> & Covnew){
    //MV t student Scale= df/(df-2) Variance
    arma::Mat<double> Scalenew=((df-2)/df)*Covnew;

    SetScale(Scalenew);
  }






  //draw new proposed value from MT(theta_t-1, df, Scale)
  arma::Col<double> DrawMVT(const arma::Col<double> & mu){

    return mu+arma::mvnrnd(Mzeros,Scale)*std::sqrt(df/arma::chi2rnd(df));

  }


  //compute the log correction factor for the proposale,
  //ie log q(theta_t|theta_prop)-q(theta_prop\theta_t)
  double logPropCorrFactor(const arma::Col<double> & loc_prop,
                           const arma::Col<double> & loc_t){
    //check symmetric PER ORA
    if(logKernel(loc_t, loc_prop)-logKernel(loc_prop, loc_t)!=0){
      throw std::runtime_error("MH not symm prop as should be");

    }

    return logKernel(loc_t, loc_prop)-logKernel(loc_prop, loc_t);


  }

  double q(const arma::Col<double> & x,
           const arma::Col<double> & mu){

    return  logconst+ logKernel(x, mu);
  }



};

class IndependentTProposal {
private:
  arma::Mat<double> Scale;    // Must be diagonal
  arma::Col<double> diagScale; // Extracted standard deviations
  double df;
  double dim;
  double logconst;  // Optional: constant part of density

public:
  IndependentTProposal(const arma::Mat<double>& scale_, double df_ = 5)
    : Scale(scale_), df(df_) {

    //  if (!Scale.is_diag()) {
    //   throw std::invalid_argument("Scale matrix must be diagonal for independent proposal.");
    //}

    diagScale = Scale.diag();  // Extract std deviations
    dim = static_cast<double>(Scale.n_cols);

    // Optional: log constant part of the t-density
    logconst = dim * (std::lgamma((df + 1.0) * 0.5) - std::lgamma(df * 0.5) - 0.5 * std::log(df * arma::datum::pi));
    logconst -= arma::sum(arma::log(diagScale));
  }

  // Draw from independent t for each dimension
  arma::Col<double> DrawMVT(const arma::Col<double>& mu) {
    arma::uword d = mu.n_elem;

    arma::Col<double> Z = arma::randn<arma::Col<double>>(d);
    arma::Col<double> U = arma::chi2rnd(df, d, 1);
    arma::Col<double> scale_factors = arma::sqrt(df / U);

    return mu + diagScale % (Z % scale_factors);
  }

  // Log kernel of sum of univariate t-densities (no constants)
  double logKernel(const arma::Col<double>& x, const arma::Col<double>& mu) {
    arma::Col<double> diff = (x - mu) / diagScale;
    arma::Col<double> squared = arma::square(diff);
    return -0.5 * arma::sum((df + 1.0) * arma::log(1.0 + squared / df));
  }

  // Symmetry check: should be zero
  double logPropCorrFactor(const arma::Col<double>& loc_prop,
                           const arma::Col<double>& loc_t) {
    double fwd = logKernel(loc_t, loc_prop);
    double rev = logKernel(loc_prop, loc_t);
    if (std::abs(fwd - rev) > 1e-10) {
      throw std::runtime_error("Independent MH proposal not symmetric as expected.");
    }
    return 0.0;
  }

  // Full log-density (if needed)
  double q(const arma::Col<double>& x, const arma::Col<double>& mu) {
    return logconst + logKernel(x, mu);
  }
};














class DiscreteXModel{
  private:// private:
    //NOTE std::vec is better than arma if we dont check bound

    arma::field<arma::Col<double>>   Adr;
    arma::field<arma::Col<double>>   logAdr;
    arma::Col<double> Adsumr;
    arma::Col<double> logAdsumr;

    std::vector<double>  loggammaAdsumr;
    arma::Col<double> sumrloggammaAdr;

    arma::Col<double>   ConstMarginalDist;


    arma::uword DD;
    //  size_t NN;
    bool is_null;
    // std::vector<arma::Col<double>>Adr;

public:
  DiscreteXModel(Rcpp::List list, arma::Col<arma::uword> MaxcatDintX, arma::uword DintX)     //Create structure: constructor

    : Adr(), logAdr(), Adsumr(), logAdsumr(), DD(0), // NN(0),
      is_null(true)  {

    if (list.containsElementNamed("a")) {
      // Initialize DD and NN based on the provided matrix in list["a"]
      arma::Mat<double> matr = as<arma::Mat<double>>(list["a"]);
      DD = matr.n_cols;
      //NON SERVE   NN = matr.n_rows;
      Adr.set_size(DD);
      logAdr.set_size(DD);
      Adsumr.set_size(DD);
      logAdsumr.set_size(DD);
      //    logprior.resize(DD);
      loggammaAdsumr.resize(DD);
      sumrloggammaAdr.set_size(DD);
      ConstMarginalDist.set_size(DD);
      is_null = false;

      for(arma::uword d=0; d<DD; d++){
        arma::Col<double> V=matr.col(d);
        arma::Col<double> logV=log(matr.col(d));

        arma::uvec indices=find(V >0);

        Adr(d)=arma::Col<double>(V.elem(indices));
        logAdr(d)=arma::Col<double>(logV.elem(indices)); //useless
        Adsumr(d)=sum(V.elem(indices));
        logAdsumr(d)=log(sum(V.elem(indices))); //useless


        //save loggamma(sum_r a_0dr) [size D]
        loggammaAdsumr.at(d)=std::lgamma(Adsumr.at(d));
        //save sum loggamma(a0dr) [size D]

        sumrloggammaAdr.at(d)=arma::sum( arma::lgamma(Adr(d)));

        //save  constant Marginal dist log( gamma(sum_r a_{0,d,r})) - sum loggamma(a0dr)
        ConstMarginalDist(d) = std::lgamma(Adsumr.at(d))-arma::sum( arma::lgamma(Adr(d)));



      }

    }else{
      throw std::runtime_error("DiscreteXModel: missing a_0,d hyperparameters");

    }




    //Check prior Discrete

    if(DD!= DintX){
      throw std::runtime_error("DiscreteXModel: Number of discrete variable DD taken from the list and those in the  Dvar main class do not match");

    }

    if(infoPrior()(0) != DintX||infoPrior().n_elem != (MaxcatDintX.n_elem+1)){


      throw std::runtime_error("DiscreteXModel: Number of discrete variable in prior not match the model' number of discrete variables");

    }else{
      for(arma::uword d=0; d<MaxcatDintX.n_elem;d++){
        if(MaxcatDintX(d)!=infoPrior()(d+1) ){
          cout << "d="<< d<< endl;
          throw std::runtime_error("DiscreteXModel: Max categories in Prior and model does not match");
        }
      }
    } ;





    //cout << "---------------------------------------------------------------------" << endl;
    //cout << "-> Initialization  DiscreteXModel class (with "<< DD<< " discrete variables)" <<endl;
    // cout << "                                                        .... OK" << endl;
    cout << "-> Initialization  DiscreteXModel class ("<< DD<< " Xd variables)...... OK"<<endl;


    //printinfo();

  }
  //#######
  //METHODS
  void printinfo(){
    for(arma::uword d=0; d<DD; d++){
      cout << "Variable d= "<< d << "with max categories "<< Adr(d).n_rows << endl;
      cout << "Probs: "<< endl;
      cout << Adr(d)  << endl;
      cout << logAdr(d)  << endl;
    }
    /*  cout << "---" << endl;

     cout << Adsumr <<   endl;
     cout << "---" << endl;
     cout << "---" << endl;

     cout << logAdsumr  << endl;*/

  }

  //method to check the coherence  between the data and the prior

  arma::Col<arma::uword> infoPrior(){
    arma::Col<arma::uword> Inf(1+Adsumr.size());
    //1st value Inf =DD number of disc variable from 1 to D: max categories
    Inf(0)=DD;

    for(arma::uword d=0; d<DD; d++){

      Inf(d+1)=Adr(d).n_rows;
    }

    return Inf;
  }

  //Methods for discrete prior predictive SINGLE Variable X_d
  double PriorPerdDisc(arma::uword const & d,arma::uword const & xdi ){
    return Adr(d)(xdi)/Adsumr(d);
  }
  //log prior predictive of in in xid
  double logPriorPerdDisc(arma::uword const & d,arma::uword const & xdi ){
    return logAdr(d)(xdi)-logAdsumr(d);
  }



  //Methods for discrete postpred predictive SINGLE variable X_d
  double logPostPerdDisc(arma::uword const & d,
                         arma::uword const & xdi,
                         double const & skdxidlessi,
                         double const & nvklessi){
    return log(Adr(d)(xdi)+skdxidlessi)-log(nvklessi+Adsumr(d));
  }

  //Marginal distribution
  double logMarginalDisc2(arma::uword const & d,
                          arma::Col<double> const &svkd,
                          // double const & svkd,
                          double const & nvk){

    double cons= std::lgamma(Adsumr(d)) - arma::sum(arma::lgamma(Adr(d)) );



    cons+= arma::sum(arma::lgamma(Adr(d)+svkd) )-std::lgamma(nvk+Adsumr(d));

    //   cout << "Inefficent complete const" << cons<<endl;

    return cons ;
  }


  //Marginal distribution
  double logMarginalDisc(arma::uword const & d,
                         arma::Col<double> const &svkd,
                         // double const & svkd,
                         double const & nvk){

    double cons=ConstMarginalDist(d);


    cons+= arma::sum(arma::lgamma(Adr(d)+svkd) )-std::lgamma(nvk+Adsumr(d));

    return cons ;
  }
  //Marginal distribution  NO const
  double logMarginalDiscNoConst(arma::uword const & d,
                                arma::Col<double> const &svkd,
                                // double const & svkd,
                                double const & nvk){

    return arma::sum(arma::lgamma(Adr(d)+svkd) )-std::lgamma(nvk+Adsumr(d)) ;
  }

};



class State  {
  private: //  private
    //Object defining the state of the chain
    arma::Mat<arma::uword> ZZ;     // Nx Nview matrix to store the view-specific allocation variables
    arma::Col<arma::uword> gammas; // vector of gammas: 1st D entries are for the disc var
    arma::Col<double>  AlphasDPs; //(Hyperparameter) DP alpha [still updated]

    
    //save initial values of gammas and cluster allocations (to be exported)
    const arma::Mat<arma::uword> InitialZZ;
    const arma::Col<arma::uword> Initialgammas;
    const arma::Col<double> InitialAlphasDPs;
    
    
    

      //DISC MODEL
      arma::Mat<arma::uword> const & Xdisc;
      arma::Col<arma::uword> const & MaxcatD;
      arma::uword const & D;

      DiscreteXModel discXmodel;
      arma::Col<arma::uword>   const  AllindicesDiscVar;





      //CONT MODEL
      //  arma::uword const & Q;
      //  arma::Mat<double> const & Xcont;

      //Y MODEL
      arma::Col<arma::uword> const & Y;
      arma::uword const & Rylessone;

      Ymodel ymodel;

      //GENERAL reference to  W*Beta transposed (Ry-1 x N)
      //constant in the state class but it can change vary in upper level
      //arma::Mat<double> const &WBETAtrans;
      //GENERAL reference to  to the vector storing the pointers of
      // W*Beta transposed (Ry-1 x N) (beta current).
      //only RefVecPtrsAllTreatWBETAtransCurrent.at(Indextreatment) need to be used in this class
      const std::vector<std::unique_ptr<arma::Mat<double> >>  &   RefVecPtrsAllTreatWBETAtransCurrent;

      //Index treatment, Used to recover   the WBETA, from the vec of pointers
      //all treatments Indextreatment: ie,
      const arma::uword   Indextreatment;

      //(structure) general
      const arma::uword Nstate;


      const arma::uword nview;
      // Non null view indeces v=1,...,Nview-1 (used for efficiency in loop)
      // arma::Col<arma::uword>   const  IndicesNonNullView;  TOGLIERE



      std::vector<std::stack<arma::uword>> ContainerNewIndices;

      arma::uword M8neal; //parameter m in Neal's 8th algorithm


      //hyper2parameters we assume same hyper-hyper parameters for all non-null view
      //for the precision parameters of DPs: alpha_v~gamma gam(shape=Hyp2parA, scale=Hyp2parB)
      double Hyp2parA, Hyp2parB;


      //Hyperparameter gamma_d  (prior P(gamma_d=v)=omega_v,d; as log )
      arma::Mat<double>  logOmega;



      //DISCRETE MODEL (structure)
      // structure to save active components (specifically it save k's,  st  n_vk>0 )
      //                 view v :     viewspecific cluster Kv : n_vk units associated to cluster K
      //NOTE the view index here has v-1
      std::vector<std::unordered_map< arma::uword, double>> ActivecompeachView;

      //structure containing the set of indices of covariates in each view
      std::vector<unordered_set<arma::uword>> Gammadisc_v;

      //structure to save the counts of for all view, cluster and discretevariables (v,k,d)
      std::unordered_map<std::tuple<arma::uword, arma::uword, arma::uword>, arma::Col<double>,
                         tuple_hash> AllS_vkd;

      //generic key, it avoid to build every time the tuple
      std::tuple<arma::uword, arma::uword, arma::uword> CentralizeKeystore;



      void ChangeCentralizeKeystore(const arma::uword& v, const arma::uword& k, const arma::uword& d) {
        std::get<0>(CentralizeKeystore) = v;
        std::get<1>(CentralizeKeystore) = k;
        std::get<2>(CentralizeKeystore) = d;

        //check
        if(AllS_vkd.count(CentralizeKeystore)==0){
          throw std::runtime_error("Key do not exist in AllSvkd");

        }
      }


      void ChangeCentralizeKeystoreNocheck(const arma::uword& v, const arma::uword& k,
                                           const arma::uword& d) {
        std::get<0>(CentralizeKeystore) = v;
        std::get<1>(CentralizeKeystore) = k;
        std::get<2>(CentralizeKeystore) = d;

      }
      //Y MODEL structure
      //structure to save the the cluster specific parameters of Y theta_k=(theta_k,1,...,theta_k,Ryless1)
      std::unordered_map<arma::uword, arma::Col<double>> Theta_ks;

      MHProposal MHthetak;

      //Acceptace unordered map k:(acccounter, total proposal) USE double since i need to divide later
      std::unordered_map<arma::uword, arma::Col<double>> AccThetak;



      //store the current Cluster specific Loglik (given current Theta_k's and betas)
      std::unordered_map<arma::uword, double> CurrentClusterLogLik;
      //store the current loglik of Y ie sum all cluster specific loglik
      double CurrentSumClusterLogLikY=0;

      //CONTINUOUS MODEL structure
      //structure containing the set of indices of covariates in each view
      std::vector<unordered_set<arma::uword>> Gammacont_v;//NON UTILIZZATA ORA







      // Posterior predictive references: Xtilde and reference to Wtildebetacurrent
      const arma::Mat<arma::uword>  & RefXdisctilde;//PostPreddisc; //predictive
      const arma::Mat<double>  & RefWtildeBetacurr; //for safety here i can use const

      const arma::uword  NPostPred;

      arma::Mat<double> ThetatildePlusWtildebetatrans;




      //sample from Gumbel
      double SampleGumbel(){
        return -log(-log(arma::randu()));
      }


      //private method dereferencing the right WBeta pointers

      const arma::Mat<double>& GetWBETATransCurr()  const {
        return *(RefVecPtrsAllTreatWBETAtransCurrent.at(Indextreatment));
      }


public:
  //CONSTRUCTOR
  State(arma::Mat<arma::uword> Zinit, //Treatment specific
        arma::Col<arma::uword> gammainit,
        arma::Col<double> initAlphaDPs,
        arma::uword m8neal, //must be >= 1
        arma::Mat<arma::uword> const &xdiscref,
        arma::Mat<double> const &xcontref,
        arma::Col<arma::uword> const & Yref,
        // arma::Mat<arma::uword> const & XPreddiscref, DA METTERE
        //Global for all treatments
        int typeofmodel,
        arma::uword  ninitview, //Treatment specific
        arma::Col<arma::uword> const & maxcardisc, // n disc variables from main class
        arma::uword const  & d, // n disc variables from main class
        arma::uword const  & ryless1, // max cat response -1
        //  arma::Mat<double> & wboldref,
        std::vector<std::unique_ptr<arma::Mat<double> >>   &   ExtRefVecPtrsAllTreatWBETAtransCurrent,
        arma::uword    IndTreat,
        Rcpp::List ListPrior,//Global for all treatments
        double LambdaTheta,
        const arma::Mat<arma::uword> & refPostPredPredictiveX,
        const arma::Mat<double> & refLinPredwbPostPred):

  ZZ(Zinit), InitialZZ(Zinit),
  AlphasDPs(initAlphaDPs), InitialAlphasDPs(initAlphaDPs),
  Nstate(Zinit.n_rows),
  gammas(gammainit), Initialgammas(gammainit),
  M8neal(m8neal),
  MaxcatD(maxcardisc),
  nview(ninitview),
  D(d),
  Xdisc(xdiscref),
  Y(Yref),
  Rylessone(ryless1),      MHthetak(arma::eye(ryless1,ryless1)*LambdaTheta ),
  ymodel(ListPrior, ryless1),
  RefXdisctilde(refPostPredPredictiveX), RefWtildeBetacurr(refLinPredwbPostPred),
  NPostPred(refPostPredPredictiveX.n_rows),
  AllindicesDiscVar(typeofmodel==1||typeofmodel==2 ?  arma::linspace<arma::Col<arma::uword>>(0, D-1, D) : 0),
  discXmodel(typeofmodel==1||typeofmodel==2 ?   ListPrior : Rcpp::List(),maxcardisc,D),
  RefVecPtrsAllTreatWBETAtransCurrent(ExtRefVecPtrsAllTreatWBETAtransCurrent), Indextreatment(IndTreat){
    //

    //cout<<IndicesNonNullView<<endl;

    //set hyper-hyepr parameters a_v, b_v (WE ASSUME a=a_v, b_v=b)
    //recall: alpha_v ~ Gamma(a_v, b_v) hyperprior
    //hyper-hyper parameters
    if (ListPrior.containsElementNamed("Hyp2AB")) {
      // Initialize DD and NN based on the provided matrix in list["a"]
      auto temphyp2parAB=as<arma::Col<double>>(ListPrior["Hyp2AB"]);
      Hyp2parA = temphyp2parAB(0);
      Hyp2parB =temphyp2parAB(1);

      /*if(temphyp2parAB.n_elem!=nview-1){
       throw std::runtime_error("hyper-hyperparameters: AlphasDPs.n_elem!=nview-1");
      } this check is wrong if we consider same heyperhyperparameters  */
    }else{
      throw std::runtime_error("hyper-hyperparameters for the hyperprior gamma(a,b) on alpha_v not given");

    }
    //NOTE  alphas DP are at begin length
    //check
    if(AlphasDPs.n_elem!=nview-1){
      throw std::runtime_error("nview init and length of initAlphaDPs did not agree");
    }




    //initialize (setting NO updating on omegas!) hyperparameters (log) gamma_d P(gamma_d=v)=omega_v,d
    if (ListPrior.containsElementNamed("HParViewallocationlog")) {
      // Initialize DD and NN based on the provided matrix in list["a"]
      logOmega = as<arma::Mat<double>>(ListPrior["HParViewallocationlog"]);

      if(logOmega.n_cols!=D||logOmega.n_rows!=nview){
        throw std::runtime_error("HParViewallocation: logOmega.n_cols!=D||logOmega.n_rows!=nview");
      }
    }else{
      throw std::runtime_error("HParViewallocation not given");

    }


    //Check Predictive X^tilde covariates DA METTERE

    //initialize ThetatildePlusWtildebeta
    ThetatildePlusWtildebetatrans.set_size(Rylessone+1, NPostPred);
    ThetatildePlusWtildebetatrans.zeros();

    //INITIALIZE DATA structures
    //INITIALIZE IF THERE ARE DISCRETE VARIABLES
    if(typeofmodel==1||typeofmodel==2){ // if we have discrete covariates
      //resize
      Gammadisc_v.resize(nview);
      ActivecompeachView.resize(nview-1);
      ContainerNewIndices.resize(nview-1);



      //CHECK ADD INIT STRUCTURES IN THE MODEL
      if(Xdisc.n_rows!=ZZ.n_rows){
        throw std::runtime_error("nrow Zinit do not match Xdisc.n_rows (= Nobs)");

      }
      if(MaxcatD.n_elem!=D){
        throw std::runtime_error("size gammainit not equal length maxcardisc( ie MaxcatD)");
      }

      if(Xdisc.n_cols!=D){
        throw std::runtime_error("length maxcardisc do not match Xdisc.n_cold  ");

      }

      if(gammas.n_elem!=D){
        throw std::runtime_error("size gamma init do not match Xdisc.n_cold  ");

      }

      if(ZZ.n_cols!=nview-1){
        throw std::runtime_error("nview init and cols of ZZ did not agree");
      }
      if(gammas.max()> nview-1)  {
        throw std::runtime_error("some indices in Gammas > than nview");
      }

      //check if Mneal 8 is >0
      if(M8neal==0 ||M8neal>1000 ){
        throw std::runtime_error("M8neal must be >0 and <1000");

      }

      //check if the max X_d > then  maximumcatDiscrete
      for(arma::uword dd=0; dd< D;dd++){
        if((Xdisc.col(dd).max()+1)> MaxcatD(dd))  {
          throw std::runtime_error("Max X_d do not match that of the maxcardisc ");
        }

      }

      //check if the number of disc pred covariates in Predictivedataset that in Xd
      if(RefXdisctilde.n_cols!=D)  {
        throw std::runtime_error("number X_dtilde covariates in predictive dataset not equals to D");
      }

      //check if the max Xtilde_d > then  maximumcatDiscrete
      for(arma::uword dd=0; dd< D;dd++){
        if((RefXdisctilde.col(dd).max()+1)> MaxcatD(dd))  {
          cout << "d" << dd <<endl;
          throw std::runtime_error("Max X_dtilde do not match that of the maxcardisc ");
        }

      }



      //check if the max Y > then  Maxcatresponse
      if(Y.max()>(Rylessone)){
        throw std::runtime_error("Max y do not match that of the Maxcatresponse ");

      }

      //check if the WBETAtrans if the correct dimension, ie Ry-1 X Nstate
      //Nothe this che only the current But since the dims of current WB is compared
      //to the WB prop this should check if there are inconsistencies
      if(GetWBETATransCurr().n_rows!=Rylessone || GetWBETATransCurr().n_cols!=Nstate){
        throw std::runtime_error("dim WBETAtrans do not match (Ryless1) or Nstate");

      }


      //initialize Gamma discrete
      for(arma::uword dd=0; dd< D;dd++){
        Gammadisc_v[gammas(dd)].insert(dd);
      }


      //initialize ActivecompeachView and ContainerNewIndices
      for(arma::uword v=0; v<ZZ.n_cols;v++){
        arma::uword maxzvplus1=1+ZZ.col(v).max();

        ContainerNewIndices.at(v).push(maxzvplus1);

        auto& MapActivCompsViewv = ActivecompeachView.at(v);
        for(arma::uword i=0; i< ZZ.n_rows;i++){
          //  if(MapActivCompsViewv.count(ZZ.at(i,v-1))==0){
          //initialize
          MapActivCompsViewv[ZZ.at(i,v)]++;
        }


      }


      // initialize unrecorded map Theta_k and CurrentClusterLogLik if Y present
      if(Y.n_elem>0){
        //Initialize Theta_k
        //

        for(auto const & k1 : ActivecompeachView.at(0)){
          //note theta_k init sample from prior
          Theta_ks.emplace(k1.first, ymodel.SampleG0thetak());

          //init Acc unordered map
          //AccThetak.emplace(k1.first,0);
          AccThetak.emplace(k1.first, arma::Col<double>(2, arma::fill::zeros));

        }

        //test for(auto const & t1 : Theta_ks){
        //  cout<<t1.first<<":"<< t1.second.t() <<endl;}

        //  printUnordlist(CurrentClusterLogLik);

        //    cout <<"PTR state ref WB"<< RefVecPtrsAllTreatWBETAtransCurrent.at(0).get() << endl;
        //   cout <<"GetWBETATransCurr \n "<< GetWBETATransCurr() << endl;
        //initialize CurrentClusterLogLik
        InitClusterLogLik(CurrentClusterLogLik,  GetWBETATransCurr()  );
      }




      //initialize CurrentSumClusterLogLikY (we can do this in various way)
      CurrentSumClusterLogLikY=ComputeLogLikYallClustersGivenBetanew(GetWBETATransCurr() );
      ////cout << "CurrentSumClusterLogLikY: " << CurrentSumClusterLogLikY <<"diff" << test1-CurrentSumClusterLogLikY<<endl;



      //Initialize All_Svkd (for all view even if a variable is not allocated to the view)
      for(arma::uword const & dd: AllindicesDiscVar){
        for(arma::uword i=0; i<Nstate; i++){
          for(arma::uword v=1; v <nview;v++){
            //   std::cout << "d" << d << "i" << i << "v"<< v <<endl;

            ChangeCentralizeKeystoreNocheck(v,ZZ(i, v-1),dd);
            //initialize new svkd if the key do not exist
            if(AllS_vkd.count(CentralizeKeystore)==0){
              AllS_vkd[CentralizeKeystore]=arma::Col<double>(MaxcatD.at(dd)).zeros();
            }

            AllS_vkd[CentralizeKeystore](xdiscref(i,dd))++;

          }

          arma::uword v=0;
          ChangeCentralizeKeystoreNocheck(v,0,dd);
          if(AllS_vkd.count(CentralizeKeystore)==0){
            AllS_vkd[CentralizeKeystore]=arma::Col<double>(MaxcatD.at(dd)).zeros();
          }

          AllS_vkd[CentralizeKeystore](xdiscref(i,dd))++;
        }
      }
      // cout << "-> Initialization  MVBPR (no global parameters) class.... OK" << endl;
      //  cout << "---------------------------------------------------------------------"<<endl;

      cout << "-> Initialization  MVBPR (no global parameters) class............ COMPLETE" <<endl;
      //cout << "-> Initialization  yModel class........................................ OK" <<endl;

    }//end init for discrete covariates

  }

  
  
  
  //METHODS
  
  //reset acceptance counter
  void ResetAccCounter(){
    for(auto & acc: AccThetak){
      acc.second.zeros();
      
    }
    //acceptanceMeanThetak=0;
  }
  
  
  
private:
  //Print map only for debug
  void printmapactivecomp(arma::uword v){
    cout<< "--"<<endl;

    // cout<< "view " << v << ":"<<endl;
    if(ActivecompeachView.at(v).size()!=0){
      for (auto s : ActivecompeachView.at(v)){
        cout<< "k_v: " << s.first  << " n_vk: "<<s.second <<endl;

      }

      cout<< "print tuple: "<< std::get<0>(CentralizeKeystore) <<
        std::get<1>(CentralizeKeystore)<<
          std::get<2>(CentralizeKeystore)<< endl;
    }

  }

  
  
  //get alphas Dp view >=1
  double GetlogAlphasDP(const arma::uword &nnview){

    return log(AlphasDPs(nnview-1));
  }

  // const double & GetlogOmega(const arma::uword &v){
  //NON USATA
  //   return logOmega(v);
  //}



  
  
  
  //Return 1st available index for a new cluster  v>=1
  arma::uword Return1stAvilLabel(const arma::uword &v){
    arma::uword  newindex= ContainerNewIndices.at(v-1).top();
    ContainerNewIndices.at(v-1).pop();

    if(!ContainerNewIndices.at(v-1).empty()){
      return newindex;
    }else{
      ContainerNewIndices.at(v-1).push(newindex+1);
      return newindex;
    }

  }


  //Remove i: this function remove the unit i from the  count in the active components
  //Specifically if n_vk=1, it remove the cluster k and return 0, if n_vk!>1,
  //n_vk is returned  and the cluster k is removed from the map, in this way we can
  //use the n_vk with k=zold, to compute posterior predictive
  //(WITHOUT the unit i-- neal 3rd algorithm)
  //Important if retured value of Removeinv_i>=2 we need to re add the zold:n_vzold
  //to the map
  double Removeinv_i(const arma::uword &v,const arma::uword &i){
    //NOTE v must be >=1
    arma::uword vminusone = v - 1;  // Calculate v-1 once

    auto& MapActivCompsViewv = ActivecompeachView.at(vminusone);
    arma::uword const &  ClusterInd = ZZ(i, vminusone);
    //                            n_vk
    if(MapActivCompsViewv.at(ClusterInd)>=2){
      //  MapActivCompsViewv[ClusterInd]--; POSSO EVITARLO LO FACCIO ALLA FINE SE znew!=zold

      double n_vzold=MapActivCompsViewv.at(ClusterInd);
      MapActivCompsViewv.erase(ClusterInd);



      return n_vzold;


    }else{//if n_vk==1 (actually here can be also 0)
      MapActivCompsViewv.erase(ClusterInd);
      //add old cluster label to the available index to use in a new cluster
      ContainerNewIndices.at(vminusone).push(ClusterInd);


      return 0;
    }
    //IMPORTANT: at this stage the obs i  is still counted in
    //in the S_vkd(x_id) for all d
  }




  //get z_v,i, for v>=1
  const arma::uword & getZvi(const arma::uword &v,const arma::uword &i){

    return ZZ(i, v-1);
  }




  //readd zodl
  void readdzoldAv(const arma::uword &v,const arma::uword &i, const double &nkold){
    //auto& MapActivCompsViewv = ActivecompeachView.at(vminusone);

    ActivecompeachView.at(v-1)[ZZ(i,v-1)]=nkold;
  }




  //Get Active components in the view v (with n_vk) DO NOT USE for v==0
  const std::unordered_map< arma::uword, double>& GetActiveCompV(const arma::uword &v ){

    return ActivecompeachView.at(v-1);
  }



  //Get indices of the discrete covariates in the view v
  const std::unordered_set<arma::uword>& GetGammaDiscV(const arma::uword &v ){

    return Gammadisc_v.at(v);
  }





  //Get success count in the K-th cluster within the v view for the d variable, for the
  //category r (NON UTILE TOGLIERE??-)
  const double & GetAllS_vkdr(const arma::uword &v, const arma::uword &k,
                              const arma::uword &d,
                              const arma::uword & r
  ){
    ChangeCentralizeKeystore(v,k,d);

    return AllS_vkd[CentralizeKeystore](r); //S_vkd,x_id
  }





  //Get success count in the K-th cluster within the v view for the d variable, for the
  //category r= xid
  const double & GetAllS_vkxid(const arma::uword &v, const arma::uword &k,
                               const arma::uword &d,
                               const arma::uword & i  ){
    ChangeCentralizeKeystore(v,k,d);

    return AllS_vkd[CentralizeKeystore](Xdisc(i,d)); //S_vkd,x_id
  }





  //Get success count vector in the K-th cluster within the v view for the d variable
  const arma::Col<double>  & GetVectS_vkd(const arma::uword &v,
                                          const arma::uword &k,
                                          const arma::uword &d ){
    ChangeCentralizeKeystore(v,k,d);

    return AllS_vkd[CentralizeKeystore]; //S_vkd,x_id
  }

  //Get success count in the K-th cluster within the v view for the d PREDICTIVE variable, for the
  //category r= xid
  const double & GetAllS_vkxidPRED(const arma::uword &v, const arma::uword &k,
                               const arma::uword xtildeid,
                               const arma::uword & d  ){
    ChangeCentralizeKeystore(v,k,d);

    return AllS_vkd[CentralizeKeystore](xtildeid); //S_vkd,x_id
  }


  

  
public:  
  //Return the current loglik of model Y (so given theta_k's and current WBeta)
  //needed to compute Beta
  double GetLogLikYallClusters(){
    
    //FOR DEBUG---- REMOVE
    double logsum=0; //debug now
    for(const auto &csl : CurrentClusterLogLik){
      logsum+=csl.second;
    }
    if(std::abs(logsum - CurrentSumClusterLogLikY) > 1e-12){
      printUnordlist(CurrentClusterLogLik);
      cout << "__________________SUM: " << logsum <<endl;


      cout << "CurrentSumClusterLogLikY: " << CurrentSumClusterLogLikY <<endl;


      throw std::runtime_error(" -----> logsum!=CurrentClusterLogLikY ");

    }//END DEBUG----

    return CurrentSumClusterLogLikY;
  }

  
  //return loglik of model Y, so given theta_k's BUT considering a new linear predictor  NEWWBETAtrans
  double ComputeLogLikYallClustersGivenBetanew(const arma::Mat<double> & NEWWBETAtrans){
   
    double LikYallClusters=0.0;
    
    //reference cluster allocation relevant view
    //this ref is okk since ZZ do not cheng and is a temp reference to col ZZ
    const auto & RefZRelevent=ZZ.col(0);
    
    //arma::mat DEBUGmatint(NEWWBETAtrans.n_cols,8);
    //DEBUGmatint.fill(-2);

    for(arma::uword i=0; i<Nstate; i++){

      LikYallClusters+=ymodel.logLikcategoricalY(Y(i), Theta_ks.at(RefZRelevent(i)),
                                                 NEWWBETAtrans.col(i));


    }



    return LikYallClusters;

  }

private:  
  

  void UpdateRecountAllocationIrrelAll(const arma::uword &v,const arma::uword &i,
                                       const arma::uword &znew){
    // size_t vminusone = v - 1;  // Calculate v-1 once
    //const arma::Row<arma::uword> & Xdisci=Xdisc.row(i); //THIS CREATE A COPY of a temp object since .row is a subview not a Row
    //const arma::Row<arma::uword>  Xdisci=Xdisc.row(i);
    
    const auto& Xdisci = Xdisc.row(i);
    
    
    
    arma::uword& zold = ZZ(i, v-1);

    //#   cout  <<endl;
    //#cout << "zold UpdateRecountAllocationIrrel: "<< zold <<endl;
    //#cout << "znew UpdateRecountAllocationIrrel: "<< znew <<endl;
    //# cout  << endl;


    auto& MapActivCompsViewv = ActivecompeachView.at(v-1);

    //is zold in the Av? 1 yes 0 no
    auto zoldinAv = MapActivCompsViewv.count(zold);


    if(zold==znew){
      //in this case zold (or znew) can be either in Av or be a new component
      //if zold==znew is in Av, do nothing; while if   zold==znew is NOT in Av
      //means the the old cluster has been destroyed but the new drown z_vi is
      //allocated to a new cluster. Hence we can recycle the old S_vkd's, we only
      //need to add the cluster in the active components
      if(zoldinAv==0){
        MapActivCompsViewv[znew]=1;
      }

    }else{//case zold!=znew
      //is znew in the Av? 1 yes 0 no
      auto znewinAv = MapActivCompsViewv.count(znew);

      if(znewinAv==1 && zoldinAv==1){
        //both znew and zold are in the active components
        //decrease nvk; k=zold
        MapActivCompsViewv[zold]--;
        //increase nvk; k=znew
        MapActivCompsViewv[znew]++;

        for(const arma::uword &d: AllindicesDiscVar){
          //decrease svkd; k=zold
          ChangeCentralizeKeystore(v,zold,d);
          AllS_vkd[CentralizeKeystore](Xdisci(d))--;
          //increase svkd; k=znew
          ChangeCentralizeKeystore(v,znew,d);
          AllS_vkd[CentralizeKeystore](Xdisci(d))++;
        }
      }else if(znewinAv==1 && zoldinAv==0){
        //znew is in the active components and zold is not  in the active components
        //n_vk already dell
        // nvk; k=znew
        MapActivCompsViewv[znew]++;
        for(const arma::uword &d: AllindicesDiscVar){

          //erase svkd ; k=zold
          ChangeCentralizeKeystore(v,zold,d);
          AllS_vkd.erase(CentralizeKeystore);
          //increase svkd; k=znew
          ChangeCentralizeKeystore(v,znew,d);
          AllS_vkd[CentralizeKeystore](Xdisci(d))++;
        }

      }else if(znewinAv==0 && zoldinAv==1){
        //znew is not the active components and zold is in the active components
        // but the old cluster was not dell

        //decrease nvk; k=zold
        MapActivCompsViewv[zold]--;
        //init new nvk
        MapActivCompsViewv[znew]=1;

        for(const arma::uword &d: AllindicesDiscVar){
          //decrease svkd; k=zold
          ChangeCentralizeKeystore(v,zold,d);
          AllS_vkd[CentralizeKeystore](Xdisci(d))--;

          //initialize svkd and nvk; k=znew
          ChangeCentralizeKeystoreNocheck(v,znew,d);

          AllS_vkd[CentralizeKeystore]=arma::Col<double>(MaxcatD(d)).zeros();
          AllS_vkd[CentralizeKeystore](Xdisci(d))++;
        }

      }else{
        cout << "Errore non dovrei mai cadere qui"<<endl;
        throw std::runtime_error("If znew and zold not in Av, znew should be reassigned to zold ");

      }

      //change the value in ZZ reassigned
      zold=znew;
    }//end if zold==znew
  }


  //Method to update all z_vi i=1,...,n with v>2
  //ONE run of Neal's 3rd algorithm (only for discrete var)
  void UPDATEIrrViewNealThirdGumbel(const arma::uword & irrv){// update Z_v,i v>=2



    //define local varibles for
    double UnlogProb;
    double  n_vkold;
    arma::uword  zold;

    //define max and argmax
    double max=-INFINITY;
    arma::uword zstar=-1;

    // Reference activecomponents
    const auto  &Activecomponent=ActivecompeachView.at(irrv-1);
    const auto & Gammadisccurrentview=Gammadisc_v.at(irrv);

    double logAlphav=GetlogAlphasDP(irrv);


    for(arma::uword i=0; i<Nstate; i++){
      //reset max and argmax
      max=-INFINITY;
      zstar=-1;


      //-----------------
      //update single z_vi
      //-----------------
      //note the k equal to zold in Activecomponent is removed
      n_vkold= Removeinv_i(irrv, i);


      //compute log denominator SI PUO EVITARE
      //  double logNless1plusAlphaDP=0;//std::log(Nless1+mcmcstate.GetAlphasDP(irrv));

      //COMPUTE UNNORMALIZED LOG PROBS FOR ALL k in A_v \ {zold}
      //for loop cluster k in active component Av less zold (within view irrv)
      for(auto const & k : Activecomponent){
        //compute the log(n_vklessi) - [- log(alpha_V+n-1) NOT REQUIRED]
        //note that Activecomponent do not contain the z_vi, n_vklessi ==n_vk
        UnlogProb=std::log(k.second);

        //Summation over post predictive discrete variables
        for(arma::uword const & d : Gammadisccurrentview){
          //here I add to UnlogProbs the log posterior pred (less i)
          //for all d in Gamma_disc_d
          //discXmodel.logPostPerdDisc= log(a_0d_xid+S_k,d,xid)-log(nk+summrad)

          UnlogProb+= discXmodel.logPostPerdDisc(d, Xdisc(i,d),
                                                 GetAllS_vkxid(irrv,k.first,d,i), k.second);


        }

        // Update max and argmax K in Activecomponent \ zold
        UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
        if (UnlogProb > max) {
          max = UnlogProb;
          zstar =k.first ;
        }

      }

      //COMPUTE UNNORMALIZED LOG PROBS FOR k=zold IF IS ACTIVE
      if(n_vkold>=2){
        //#cout<<"HEREE" << i <<endl;

        double nvkless1 = n_vkold-1;
        //compute the log(n_vklessi) [- log(alpha_V+n-1) NOT REQUIRED]here n_vklessi!=n_vk
        UnlogProb=std::log(nvkless1);
        zold=ZZ(i, irrv-1) ;

        //Summation over post predictive discrete variables
        for(arma::uword const & d : Gammadisccurrentview){
          // add log post predictive (less i) for k=zold if in A_v
          //discXmodel.logPostPerdDisc= log(a_0d_xid+(S_k,d,xid - 1))-log(nk-1+summrad)

          UnlogProb+= discXmodel.logPostPerdDisc(d, Xdisc(i,d),
                                                 GetAllS_vkxid(irrv, zold ,d,i)-1, nvkless1);

        }
        //IMPORTANT STEP: re insert zold cluster in A_v if zold is active
        //namely n_vkold>=2
        readdzoldAv(irrv,i,n_vkold);

        // Update max and argmax for K=zold if k is ACTIVE
        UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
        if (UnlogProb > max) {
          max = UnlogProb;
          zstar =zold ;
        }
      }


      //COMPUTE UNNORMALIZED LOG PROBS FOR NEW CLUSTER
      // compute log(alpha_v)  [- log(alpha_V+n-1) NOT REQUIRED]
      UnlogProb=logAlphav;

      //Summation over prior predictive discrete variables
      for(arma::uword const & d : Gammadisccurrentview){
        // add log prior predictive (less i)
        UnlogProb+= discXmodel.logPriorPerdDisc(d, Xdisc(i,d));

      }

      // Update max and argmax for K=znew cluster
      UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
      if (UnlogProb > max) {
        max = UnlogProb;
        //since K=znew cluster is the last unnormalizedprobs so if
        //here UnlogProb > max means that we drawn znew, so we give it the
        //1st available index
        zstar =Return1stAvilLabel(irrv) ;
      }


      //UPDATE Z_vi. Note that internally manage the computetion of the update state, eg ZZ(i, v-1) 
      UpdateRecountAllocationIrrelAll(irrv,i,zstar);

    }


  }




  //Method to update all z_vi i=1,...,n with v>2
  //ONE run of Neal's 8rd algorithm (only for discrete var)
  void UPDATERelViewNealEightGumbel(){// update Z_v,i v>=2
    // set view 1
    const arma::uword & relv=1;

    // de referencing  WB at the begin of this method
    arma::Mat<double> const & WBetacurr =GetWBETATransCurr();





    //Neal 8 specificity
    //local unorderedmap save loglik clusterspecific BETA theta old (z_new)
    //  std::unordered_map<arma::uword, double> ClusterLogLikBetaThetaoldZnew;
    arma::Mat<double> Mmatrix;
    CurrentClusterLogLik.clear();
    CurrentSumClusterLogLikY=0;
    //  std::unordered_map<arma::uword, double>().swap(ClusterLogLikBetaThetaoldZnew);


    //define local varibles for gumbel and to store p(y_i|theta_zstar)
    double UnlogProb;
    double  n_vkold;
    // arma::uword  zold; not here defined leter
    double  LoglikYigivenZstar;


    //define max and argmax
    double max=-INFINITY;
    arma::uword zstar=-1;
    arma::uword IndexAvZoldNew; // will indicate if the zstar is
    //one of k in Av (if 0), if ==1 or 2  zstar= 1s available index for a new grup
    // if 1 theta_zstar=old theta if 2  theta_zstar=one of the m (or m-1) aux parameters

    arma::uword mstar; //save index of the sampled auxparameter



    // Reference activecomponents view 1
    const auto  &Activecomponent=ActivecompeachView.at(relv-1);
    const auto & Gammadisccurrentview=Gammadisc_v.at(relv);

    //save DP_ hyperparameter alpha_1 view1
    double logAlphav=GetlogAlphasDP(relv);
    //->*/
    //Reference to Zrelevant view
    //arma::Col<arma::uword> const &Zrel=ZZ.col(relv-1); //CREATE A REF TO A TEMP (copyed) OBJ
    auto const &Zrel=ZZ.col(relv-1);



    //-----------------
    //UPDATE ALL Z_1,i
    //-----------------
    //-> /*
    for(arma::uword i=0; i<Nstate; i++){
      //reset max and argmax
      max=-INFINITY;
      zstar=-1;
      IndexAvZoldNew=0;
      mstar=-1; //actually we can skip this since IndexAvZoldNew!=2 i do not need this
      LoglikYigivenZstar=0;

      //-----------------
      //update single z_vi
      //-----------------
      //note the k equal to zold in Activecomponent is removed
      //note that theta_k with k=zold still present (ie not removed) from theta_ks
      n_vkold= Removeinv_i(relv, i);

      //current value of of z_vi create a copy
      arma::uword const  zold=Zrel(i);


      //save ref to yi and w_i'beta (is a vector becuse we have Ry-1 Linpred)
      //arma::Col<double> const& WitBeta=WBetacurr.col(i); THIS CREATE A TEMP COPY
      
      auto const& WitBeta=WBetacurr.col(i);
      
      
      arma::uword const & yi=Y(i);


      //compute log denominator SI PUO EVITARE
      //  double logNless1plusAlphaDP=0;//std::log(Nless1+mcmcstate.GetAlphasDP(irrv));

      //COMPUTE UNNORMALIZED LOG PROBS FOR ALL k in A_v \ {zold}
      //for loop cluster k in active component Av less zold (within view irrv)
      for(auto const & k : Activecomponent){
        //compute the log(n_vklessi) + logP(y_i|theta_k w_i%*%beta),  [- log(alpha_V+n-1) NOT REQUIRED]
        //note that Activecomponent do not contain the z_vi, thefore here n_vklessi ==n_vk
        auto LogLikYigivenK=ymodel.logLikcategoricalY(yi,
                                                      Theta_ks.at(k.first), WitBeta);

        UnlogProb=std::log(k.second)+LogLikYigivenK;

        //Summation over post predictive discrete variables
        for(arma::uword const & d : Gammadisccurrentview){
          //here I add to UnlogProbs the log posterior pred (less i)
          //for all d in Gamma_disc_d
          //discXmodel.logPostPerdDisc= log(a_0d_xid+S_k,d,xid)-log(nk+summrad)
          UnlogProb+= discXmodel.logPostPerdDisc(d, Xdisc(i,d),
                                                 GetAllS_vkxid(relv,k.first,d,i), k.second);


        }

        // Update max and argmax K in Activecomponent \ zold
        UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
        if (UnlogProb > max) {
          max = UnlogProb;
          zstar =k.first ;
          //save p(y_i|theta_zstar, b)
          //note if we save the sum over all i, in the groups we can avoid to
          //recompure the clusterspecific loglik
          LoglikYigivenZstar=LogLikYigivenK;
        }

      }
      //COMPUTE UNNORMALIZED LOG PROBS FOR k=zold IF IS ACTIVE
      if(n_vkold>=2){
        //specificity neal's 8th algorithm
        //SAMPLE new M parameters theta_k from G_0y to be used leter
        Mmatrix=ymodel.SampleG0thetak(M8neal);

        double nvkless1 = n_vkold-1;

        //  zold=ZZ(i, relv-1) ; //current value of of z_vi

        //compute the log(n_vklessi)+logP(y_i|theta_k w_i%*%beta) [- log(alpha_V+n-1) NOT REQUIRED]
        //here n_vklessi!=n_vk
        auto LogLikYigivenK=ymodel.logLikcategoricalY(yi,
                                                      Theta_ks.at(zold), WitBeta);

        UnlogProb=std::log(nvkless1)+LogLikYigivenK;

        //Summation over post predictive discrete variables
        for(arma::uword const & d : Gammadisccurrentview){
          // add log post predictive (less i) for k=zold if in A_v
          //discXmodel.logPostPerdDisc= log(a_0d_xid+(S_k,d,xid - 1))-log(nk-1+summrad)

          UnlogProb+= discXmodel.logPostPerdDisc(d, Xdisc(i,d),
                                                 GetAllS_vkxid(relv, zold ,d,i)-1,
                                                 nvkless1);

        }
        //IMPORTANT STEP: re insert zold cluster in A_v if zold is active
        //namely n_vkold>=2
        readdzoldAv(relv,i,n_vkold);

        // Update max and argmax for K=zold if k is ACTIVE
        UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
        if (UnlogProb > max) {
          max = UnlogProb;
          zstar =zold ;
          //save p(y_i|theta_zold, b)
          //note if we save the sum over all i, in the groups we can avoid to
          //recompure the clusterspecific loglik
          LoglikYigivenZstar=LogLikYigivenK;

        }
      }else{
        //COMPUTE UNNORMALIZED LOG PROBS FOR k=zold IF IS NOT ACTIVE
        //if zold is a singleton cluster the 1st prob is computed as
        //log p(y1|theta_zold, beta)
        if(n_vkold==0){ //eventualmente si può togliere
          //SAMPLE new M parameters theta_k from G_0y to be used later
          Mmatrix=ymodel.SampleG0thetak(M8neal-1);


          //log p(y1|theta_zold, beta)
          auto LogLikYigivenK=ymodel.logLikcategoricalY(yi,
                                                        Theta_ks.at(zold), WitBeta);


          // compute log(alpha_v)-log(M)+ log p(y1|theta_zold, beta) [- log(alpha_V+n-1) NOT REQUIRED]
          UnlogProb=logAlphav-log(M8neal)+LogLikYigivenK;

          for(arma::uword const & d : Gammadisccurrentview){
            // add log prior predictive (less i)
            UnlogProb+= discXmodel.logPriorPerdDisc(d, Xdisc(i,d));

          }


          // Update max and argmax for K=zold if k is ACTIVE
          UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
          if (UnlogProb > max) {
            max = UnlogProb;
            //if sampled this zstar MUST BE zold but for all new cluster we need
            //sampled also all the other m -1 parameters

            IndexAvZoldNew=1;

            //save p(y_i|theta_zold, b)
            //note if we save the sum over all i, in the groups we can avoid to
            //recompute the clusterspecific loglik
            LoglikYigivenZstar=LogLikYigivenK;

          }

        }else{
          throw std::runtime_error("Non devo mai cadere qui: neal8 n_vkold<2 and n_vkold!=0");
          //NOTE that Removeinv_i returns n_zoldi if n_zoldi>=2 but if the old cluster is a singleton 
          //ie n_zoldi=1 then Removeinv_i=0 not 1   
        }

      }


      //COMPUTE UNNORMALIZED LOG PROBS FOR NEW CLUSTER AUXILIARY
      //PARAMETERS M8neal or M8neal-1
      for(arma::uword mindex = 0; mindex < Mmatrix.n_cols; mindex++){


        //log p(y1|theta_m~G0, beta)
        auto LogLikYigivenK=ymodel.logLikcategoricalY(yi,
                                                      Mmatrix.col(mindex), WitBeta);


        // compute log(alpha_v)-log(M)+log p(y1|theta_zold, beta)  [- log(alpha_V+n-1) NOT REQUIRED]
        UnlogProb=logAlphav-log(M8neal)+LogLikYigivenK;


        //Summation over prior predictive discrete variables
        for(arma::uword const & d : Gammadisccurrentview){
          // add log prior predictive (less i)
          UnlogProb+= discXmodel.logPriorPerdDisc(d, Xdisc(i,d));

        }

        // Update max and argmax for K=znew cluster
        UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
        if (UnlogProb > max) {
          max = UnlogProb;
          IndexAvZoldNew=2;
          mstar = mindex;
          //save p(y_i|theta_zold, b)
          //note if we save the sum over all i, in the groups we can avoid to
          //recompute the clusterspecific loglik
          LoglikYigivenZstar=LogLikYigivenK;
        }



      }




      //Generate zstar if sampled one auxiliary parameter and
      //rearrange the state for thetak
      if(IndexAvZoldNew!=0){
        //    cout << "----------->>>NEW THETA SAMPLED<<<-"<<endl;
        //Sample zstar for new cluster
        zstar=Return1stAvilLabel(relv);

        if(n_vkold>=2){//old cluster not destroyed
          //add m-th aux parameter to the state
          Theta_ks.emplace(zstar, Mmatrix.col(mstar));
          if(IndexAvZoldNew==1){throw std::runtime_error("!!!IndexAvZoldNew==1"); }


          // AccThetak.emplace(zstar,0);
          AccThetak.emplace(zstar, arma::Col<double>(2, arma::fill::zeros));




        }else{//old cluster destroyed
          if(IndexAvZoldNew==2){
            //  cout << "----------->>>An old THETA replaced<<<-"<<endl;

            if(zstar!=zold){throw std::runtime_error("!!!flagSampleZold"); }

            //replace theta old with theta m aux
            Theta_ks.at(zstar)=Mmatrix.col(mstar);
            //can i make more efficent Theta_ks[zstar] ? Yes


            //reset theta_k in Acc
            AccThetak.at(zstar).zeros();
          }

        }

      }else{
        if(n_vkold==0){
          //   cout << "----------->>>An old THETA erased<<<-"<<endl;

          //need to erase theta old
          Theta_ks.erase(zold);

          //erase also the theta_k from Acc
          AccThetak.erase(zold);

        }
      }


      //save clusterspecificlik: Sum_i in C_1k P(y_i| theta_k^[s], beta^[s-1])
      CurrentClusterLogLik[zstar]+=LoglikYigivenZstar;





      //this method recount all Svkd and update the A_v  and change Zvi
      //Must be at the end of the for i
      UpdateRecountAllocationIrrelAll(relv,i,zstar);


    }//end Z1 updating

    //



    //-------------------
    //UPDATE ALL theta_k
    //-------------------
    //def variable to save the loglikclusterk given the proposed theta
    double LoglikYigivenThetaProp;




    //perhaps not so efficient :(
    for(auto  & kthetak : Theta_ks){// Theta_ks contains all current thetak

      const arma::uword k = kthetak.first;
      const arma::vec& ThetakCurrent = kthetak.second;
      //   cout <<"proposed::k " << k<< endl;

      //Draw from the proposal q(.| theta_t)
      auto Thetakprop=MHthetak.DrawMVT(ThetakCurrent);

      //compute loglik cluster k for theta_k=theta_k proposed
      LoglikYigivenThetaProp=0.0;

      
      //maybe i can create a map (need to be recomputed ) outside for kthetak and a for in in this map[k]
      
      for(arma::uword i=0; i<Nstate; i++){
        //    cout <<"proposed::i " << i << " zi" << Zrel(i)<<endl;
        if(Zrel(i)==k){
          //  cout <<"Zrel(i)==k " <<endl;
          LoglikYigivenThetaProp+=ymodel.logLikcategoricalY(Y(i),
                                                            Thetakprop,
                                                            WBetacurr.col(i));

 

        }
      }




      //MH step
      double logr=LoglikYigivenThetaProp+ymodel.logPriorindtstud(Thetakprop)- //logpost prop
        CurrentClusterLogLik.at(k)-ymodel.logPriorindtstud(ThetakCurrent)+//logpost current+
        MHthetak.logPropCorrFactor(Thetakprop,ThetakCurrent); //here should be 0

      if(log(arma::randu())<logr){//accept proposed with prob min(1, exp(logr))
        //cout << "accepted theta: "<<  ThetakCurrent.t() << "newloglik k"<< test[k]<<endl;
        //cout << "-->ACC thetak---" <<endl;
        // std::cout << "acce"<<accmean<<endl;

        //update the Theta_k in the unorderedmap Theta_ks
        kthetak.second=Thetakprop;

        //update the acceptance
        //acceptanceMeanThetak+=accmean;
        AccThetak.at(k)++;


        //update the k-th cluster specific likelihood, now will be that computed with theta_prop
        //i.e., LoglikYigivenThetaProp
        CurrentClusterLogLik.at(k)=LoglikYigivenThetaProp; //(QUESTO MAGARI NON SERVE PIU)

        CurrentSumClusterLogLikY+=LoglikYigivenThetaProp;
      }else{
        AccThetak.at(k).at(1)++;
        
        CurrentSumClusterLogLikY+=CurrentClusterLogLik.at(k);

      }





    }

  }






  //Method to update all gamma_d for d=1,...,D (only for discrete variables)
  void UPDATEGammasDisc(){// update Z_v,i v>=2

    double UnlogProb;

    //define max and argmax  (gammaold)
    double max;
    arma::uword gammastar;

    for(const arma::uword & d: AllindicesDiscVar ){
      arma::uword & gammaold=gammas(d);

      //reset max and argmax
      max=-INFINITY;
      gammastar=-1;

      //COMPUTE for gamma_d, probs to be equal to v (1,...L-1) non-null view
      for(arma::uword view=1; view<nview; view++ ){

        // log prior gamma_d (SAME FOR ALL d)
        UnlogProb= logOmega(view,d);


        for(auto const & k : GetActiveCompV(view)){

          //add marginal lik x_id: lgamma(sumr a_0dr) - sum_r lgamma(a_0dr) + sim_r lgamma(a_0dr+s_vkdr)-lgamma(n_vk+sumr a_odr)
          UnlogProb+=discXmodel.logMarginalDisc(d,GetVectS_vkd(view, k.first, d), k.second);

        }


        // Update max and argmax for view =1,... L-1 (non null)
        UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
        if (UnlogProb > max) {
          max = UnlogProb;
          gammastar =view ;
        }

      }

      //COMPUTE for gamma_d, probs to be equal to v =0  null view

      UnlogProb= logOmega(0,d)+discXmodel.logMarginalDisc(d,GetVectS_vkd(0, 0, d), Nstate);
      // cout<<  logOmega(0,d) ;
      //NOTE THAT THIS DONT CHANGE WITh ITErATIONS we MAY PRE COMPUTE


      UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
      if (UnlogProb > max) {
        max = UnlogProb;
        gammastar =0 ;
      }


      //Update the state
      if(gammastar!=gammaold){
        Gammadisc_v.at(gammaold).erase(d);
        Gammadisc_v.at(gammastar).emplace(d);
        gammaold=gammastar;
      }

    }

  }


  //method to update alpha parameters of DPs
  void UPDATEalphaDP(const arma::uword  & v){
    //Escobar and West, 1995
    double u, pi_Xi, aplusKlessone, blesslogXi, aux1, aux2;
    arma::uword vless1;

    vless1=v-1;
    //draw Xi ~beta(alpha_vOLD, n)
    //arma do not have beta: we use 2 gamma
    //note that randg(shape, scale) while rgamma(shape, rate=1/scale)
    aux1=arma::randg<double>(arma::distr_param(AlphasDPs(vless1)+1.0, 1.0));
    aux2=arma::randg<double>(arma::distr_param(Nstate, 1.0));

    //compute b-log(Xi)
    blesslogXi=Hyp2parB-log(aux1/(aux1+aux2));

    //compute a0+K-1
    aplusKlessone=Hyp2parA-1+static_cast<double>(ActivecompeachView.at(vless1).size());
    // cout << ActivecompeachView.at(v-1).size() <<endl;

    //compute pi_xi
    pi_Xi=aplusKlessone/(aplusKlessone+Nstate*blesslogXi);

    //sample c=ber(pi_xi)
    u = arma::randu();
    if(u < pi_Xi){//if c=1
      // UPDATING gamma(a+k=aplusKlessone+1, 1/blesslogXi)
      AlphasDPs(vless1)=arma::randg<double>(arma::distr_param(aplusKlessone+1, 1/blesslogXi));;
    }else{//if c=0
      // UPDATING gamma(a+k=aplusKlessone, 1/blesslogXi)
      AlphasDPs(vless1)=arma::randg<double>(arma::distr_param(aplusKlessone, 1/blesslogXi));;

    }
  }

  
  //These methods must be accessible outside the functions
public:
  //get accepted  AccThetak
  const  std::unordered_map<arma::uword, arma::Col<double>> & getAccThetak(){
    
    return AccThetak;
  }

  //method to perform 1 run of MCMC sample (with state updating)
  void SAMPLEROneRun(){

    //UPDATE Relevant view=1
    UPDATERelViewNealEightGumbel();
    UPDATEalphaDP(1);



    //UPDATE Irrelevan view=2,...,Nview-1
    for(arma::uword irrelevantv=2;irrelevantv<nview;  irrelevantv++){
      UPDATEIrrViewNealThirdGumbel(irrelevantv);
      UPDATEalphaDP(irrelevantv);
    }
    //UPDATE gammas
     UPDATEGammasDisc(); //RIMETTERE SE TOLTO DALLA FUNZIONE

  }



  //method to Compute posterior predictive distribution response
  //arma::Col<arma::uword>
  inline void COMPUTEPostPredResponse(arma::umat   & YPostpred, arma::uword mcmcind ){

    //relevant view index
    arma::uword  relv=1;

    //index sample ztilde
    double max=-INFINITY;
    arma::uword ztildestar=-1;
    double UnlogProb=0;


    //references Active clusters in relevant view and Gammas
    const auto  &ActiveComponentRelevantView=ActivecompeachView.at(relv-1);
    const auto & GammadiscRelevantView=GetGammaDiscV(relv);

    const double logAlphav=GetlogAlphasDP(relv);

    const auto & Thetarelevant=Theta_ks;

    //Note GetGammaDiscV(relv)== Gammadisc_v.at(relv)
    //Compute W*Beta



    //Compute Ztilde and populate the matrix theta_tilde for each units
    for(arma::uword itilde=0; itilde<NPostPred; itilde++){

      //reset max and argmax
      max=-INFINITY;
      ztildestar=-1;



      //COMPUTE UNNORMALIZED LOG PROBS FOR ALL ACTIVE k in A_v
      for(auto const & k : ActiveComponentRelevantView){
        //compute the log(n_k) - [- log(alpha_V+n) NOT REQUIRED]
        //note that k.second = n_k instead k.first is the index
        UnlogProb=std::log(k.second);

        //Summation over post predictive discrete variables
        for(arma::uword const & d : GammadiscRelevantView){
          //here I add to UnlogProbs the log posterior
          //for all d in Gamma_disc_d
          //discXmodel.logPostPerdDisc= log(a_0d_xid+S_k,d,xid)-log(nk+summrad)

          UnlogProb+= discXmodel.logPostPerdDisc(d, RefXdisctilde(itilde,d),
                                                 GetAllS_vkxidPRED(relv,k.first,RefXdisctilde(itilde,d),d), k.second);


        }

        // Update max and argmax Ktilde in Activecomponent
        UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise
        if (UnlogProb > max) {
          max = UnlogProb;
          ztildestar =k.first ;
        }
      }

      //COMPUTE UNNORMALIZED LOG PROBS FOR NEW CLUSTER
      // compute log(alpha_v)  [- log(alpha_V+n) NOT REQUIRED]
      UnlogProb=logAlphav;

      //Summation over prior predictive discrete variables
      for(arma::uword const & d : GammadiscRelevantView){
        // add log prior predictive (less i)
        UnlogProb+= discXmodel.logPriorPerdDisc(d, RefXdisctilde(itilde,d));

      }
      // Update max and argmax for Ktilde=znew cluster
      UnlogProb+=SampleGumbel(); // Perturbation UnlogProb std gumbel noise


      //IF UnlogProb > max means that the sampled Ztilde is a new cluster,
      //otherwise is one in A_v, specifically the ztildestar.
      //So we need to populate the matrix Thetatilde accordingly

      if (UnlogProb > max) {//case: Ktilde=znew cluster
        //  max = UnlogProb;
        //since K=znew cluster is the last unnormalizedprobs we can save a
        //draw from G0_theta in the matrix thetatilde

        ThetatildePlusWtildebetatrans(arma::span(0, Rylessone-1), itilde)=ymodel.SampleG0thetak();


        //ThetatildePlusWtildebeta(itilde, 0)=1000000;
        //Ztilde(itilde, mcmcind)=1000000;

      }else{//case: Ktilde on in A_v,
        //cout << "B" << endl;
        ThetatildePlusWtildebetatrans( arma::span(0, Rylessone-1), itilde)=Thetarelevant.at(ztildestar);


        //Ztilde(itilde, mcmcind)=ztildestar;
      }
    }
    //IMPORTANT setp:  reset to zero the last row, the row of the reference category of the response
    //because ThetatildePlusWtildebetatrans is a global object of this class and the function
    //SampleResponse take the reference of a matrix, add gumbel noise and return the sampled response
    //thefrore it modify the original matrix 
    ThetatildePlusWtildebetatrans.row(Rylessone).zeros();




    //here ThetatildePlusWtildebeta have all theta_ztilde, we need to add Wtildbeta
    //which is common across treatmens and NEED TO BE updated to the current beta
    //BEFORE calling this function with the method in the GlobalParameters class
    ThetatildePlusWtildebetatrans.rows(0, Rylessone-1) += RefWtildeBetacurr;



    //Sample the response and return
    YPostpred.col(mcmcind)=ymodel.SampleResponse(ThetatildePlusWtildebetatrans);

  }




  //method to save all parameters and compute posterior predictive
  void SaveParamsandPostPred(arma::ucube & Zsave,
                             arma::umat & GammasSave,
                             arma::cube & Thetaus,
                             arma::mat & Alphasave,
                             arma::umat   & YPostpred,
                             arma::uword mcmcind){
    
    Zsave.slice(mcmcind)=ZZ;
    
    GammasSave.col(mcmcind)=gammas;
    
    
    Alphasave.col(mcmcind)=AlphasDPs;
    
    
    COMPUTEPostPredResponse(YPostpred, mcmcind);
    
    
    for(arma::uword  i=0; i<Nstate; i++){
      Thetaus.slice(mcmcind).col(i)=Theta_ks[getZvi(1,i)];
      
    }
    
  }
  
  
   
  
  
private:

  //Methods to save results
  void SaveThetaunitspecific(arma::cube  & Thetaus, arma::uword mcmcind){

    for(arma::uword  i=0; i<Nstate; i++){
      Thetaus.slice(mcmcind).col(i)=Theta_ks[getZvi(1,i)];

    }
  }








  void SaveAlphaDPs(arma::mat & Alphasave , arma::uword mcmcind){
    Alphasave.col(mcmcind)=AlphasDPs;
  }

  void SaveViewallocantos(arma::umat & GammasSave, arma::uword mcmcind){
    GammasSave.col(mcmcind)=gammas;
  }
  void SaveClusterallocantos(arma::ucube & Zsave, arma::uword mcmcind){
    Zsave.slice(mcmcind)=ZZ;

  }


  ///checker State view1
  //after Update recount
  void checkerState1(){
    //check k:nk
    std::unordered_map<arma::uword, double> auxAv;
    for(auto & zi : ZZ.col(0)){
      auxAv[zi]+=1;
    }

    //check if auxAv and Av
    auto const Avstate =GetActiveCompV(1);
    if(auxAv.size()!=Avstate.size()){
      throw std::runtime_error("auxAv.size()!=GetActiveCompV(1).size()");

    }else{

      for(const auto& k : auxAv){
        if(Avstate.count(k.first)==0){
          throw std::runtime_error("k of ZZ not match k in state Av");
        }

        if(Avstate.at(k.first)!=k.second){
          throw std::runtime_error("n_k of ZZ not match n_k in state Av");
        }
      }

    }

    //check if clusters label in Av match those in Thetak
    if(Avstate.size()!=Theta_ks.size()){
      throw std::runtime_error("Avstate.size()!=Theta_ks.size()");
    }else{
      for(const auto& k : Avstate){
        if(Theta_ks.count(k.first)==0){
          throw std::runtime_error("Theta_ks.count(k.first)==0");
        }
      }

    }
    //check if clusters label in Av match those in Thetak
    if(Avstate.size()!=Theta_ks.size()){
      throw std::runtime_error("Avstate.size()!=Theta_ks.size()");
    }else{
      for(const auto& k : Avstate){
        if(Theta_ks.count(k.first)==0){
          throw std::runtime_error("Theta_ks.count(k.first)==0");
        }
      }

    }



  }

  //auxiliary method to compute the log lik cluster used in constructor
  void InitClusterLogLik(std::unordered_map<arma::uword, double> &UnMap,
                         const arma::Mat<double> & wbt){

    arma::uword v=1;
    //check if the UnMap is empty otherwise
    if(!UnMap.empty()){
      cout<< "Erese Unordered map" << endl;
      UnMap.clear();
    }


    //compute the loglik for all  clusters
    const auto Z=ZZ.col(v-1);


    for(arma::uword i=0; i<Nstate; i++){
      UnMap[Z(i)]+=ymodel.logLikcategoricalY(Y(i),
              Theta_ks.at(Z(i)), wbt.col(i));


    }
  }



  //checker loglik cluster (after Z1 updating)
  void checkerclusterloklik(std::unordered_map<arma::uword, double> CLLBTON){
    //
    auto const Avstate =GetActiveCompV(1);

    //check if cluster log have same keys of Av
    if(Avstate.size()!=CLLBTON.size()){
      throw std::runtime_error("Avstate.size()!=ClusterLogLikBetaThetaoldZnew.size()");
    }else{
      for(const auto& k : Avstate){
        if(CLLBTON.count(k.first)==0){
          throw std::runtime_error("CLLBTON.count(k.first)==0");
        }
      }
    }

    auto const wbt=GetWBETATransCurr();

    for(const auto& k : Theta_ks  ){
      double naivelikK=0;
      for(arma::uword ii=0; ii<Nstate; ii++){
        if(ZZ(ii,0)==k.first){
          naivelikK+=
            ymodel.logLikcategoricalY(Y(ii),
                                      Theta_ks.at(k.first), wbt.col(ii));


        }
      }

      if(naivelikK!=CLLBTON.at(k.first)){
        cout<< naivelikK-CLLBTON.at(k.first)<<endl;
        throw std::runtime_error("naivelikK!=CLLBTON(k.first)");




      }


    }
  }


  //// Checker All S_vkd
  void checkerAllsvkd(){
    cout << "+++view 0+++"<<endl;

    //view 0 sod
    for(arma::uword d=0;d<D; d++){
      ChangeCentralizeKeystore(0,0,d);

      for (arma::uword val = 0; val < MaxcatD(d); ++val) {

        arma::uvec indices = arma::find(Xdisc.col(d) == val);

        auto s0d=indices.n_elem;

        if(AllS_vkd[CentralizeKeystore](val)!=s0d){
          throw std::runtime_error("errore view 0");
        }

      }

      if(arma::sum(AllS_vkd[CentralizeKeystore])!=Xdisc.n_rows){
        throw std::runtime_error("errore view 0 not sum to n");


      }


    }




    for(arma::uword v=1; v<nview; v++){
      cout << "+++view z+++"<<endl;

      unordered_map<arma::uword,double> Klust;

      for(auto & zi : ZZ.col(v-1)){
        Klust[zi]++;


      }

      //check if active component are choearent

      for(auto & k : Klust){
        if(ActivecompeachView.at(v-1).count(k.first)==0){
          throw std::runtime_error("errore 1");
        }else{
          if(k.second!=ActivecompeachView.at(v-1)[k.first]){

            throw std::runtime_error("errore 1n");
          }

        }
      }
      for(auto & k :ActivecompeachView.at(v-1) ){
        if(Klust.count(k.first)==0){
          throw std::runtime_error("errore 2");
        }else{
          if(k.second!=Klust[k.first]){


            throw std::runtime_error("errore 2n");
          }

        }
      }




      cout << "+++view view v+++"<<endl;

      //skvd
      for(arma::uword d=0;d<D; d++){
        arma::uvec dind(1);

        dind(0)=d;

        for(auto & k :ActivecompeachView.at(v-1) ){
          ChangeCentralizeKeystore(v,k.first,d);

          arma::uvec indicesZk = arma::find(ZZ.col(v-1) == k.first);


          // Extract subview



          for (arma::uword val = 0; val < MaxcatD(d); val++) {

            arma::uvec indices = arma::find(  Xdisc.submat(indicesZk,dind)== val);

            auto svkdval=indices.n_elem;

            if(AllS_vkd[CentralizeKeystore](val)!=svkdval){
              cout << "xid"<< val<<endl;
              cout << "svkdval"<< svkdval<<endl;
              cout << "AllS_vkd"<<   AllS_vkd[CentralizeKeystore](val)<<endl;


              throw std::runtime_error("errore view view v k");
            }

          }
          if(arma::sum(AllS_vkd[CentralizeKeystore])!=k.second){
            throw std::runtime_error("errore view 0 not sum to n");


          }

        }


      }


    }
    cout << "Svkdr OK!"<<endl;

  }


};






class GlobalParameters{
private:
  //SOME REFERENCES
  // some reference to ``global quantities''
  arma::uword const & Rylessone; //reference max number of categories  response model
  //arma::uword RylessTwo, RylessOne;

  arma::uword const & ntreatments; //Number of treatments 0... T-1

  arma::uword const nProgn; //Number prognostic

  //HYPERPARAMETERS
  //Hyperparameters for beta, NOTE we have  (Ry-1)xP parameters, which are modeled as indep t student
  //conveniently, we store the (Ry-1)xP parameters beta in a matrix:
  // if W is the n x P matrix containinng the prognostic variables, and the matrix
  // Bis a p x Ry-1 matrix, so the r column represent the effect on the
  //log odds of category r vs Ry of the p prognostic variable, then the linear
  //predictor regarding the prognostic variables, for all units is given by W*B,
  //However is is more convinient to work with its transposed (B'*W') [Ry-1 x n]
  //this because when we work with the unit i, we need w_i'*B_1,...,w_i'*B_Ry-1,
  //and taking the i-th column of  (B'*W') it is more efficent and return a col vec.
  //thefore it core convenient to consider B^~=B' [Ry-1 x P]
  //Hyperparameters for beta
  arma::Mat<double> beta0, scale, df;







  //auxiliary data structure to efficient computation of the prior quantities (logprior)
  arma::Mat<double> inv_dfscale2,  NegdfPlusOnedivTwo;
  double logNorConst;


  //store as reference all Prognostic variables for all treatments the
  //arma::Cube<double> const & Walltreatments;
  const std::vector<std::unique_ptr<const arma::Mat<double> >>  & WtransAllTreatments_ptr;

  //store as reference State of the MVBPR por each treatments
  const std::vector<std::unique_ptr<State>>  & STATEalltreatmets_ptr;

  //current Beta and proposed beta mat
  arma::Mat<double> BetasMatCurr;
  arma::Mat<double> BetasMatProp;


  //matrix to save online sum_s Beta_r, Beta_r^t and sum_s Beta_r
  std::optional< arma::Cube<double> > CubeAllcarAccumBetarBetartrans;
  std::optional< arma::Mat<double> > MatrixAllCatAccumBetar;
  double SampleSizeVariance=0;


  //store the curren loglikY associated to the current BETA and store the
  //log prior associated to the current BETA
  // double logpriorcurrent; [STO USANDO UNA VARIABILE LOCALE IN UPDATE--SE VOLESSI RISPARMIARE DEVO USARE QUESTA GLOBALE]


  //Matrices [Ry-1 x Nstate]proposed (WBetaprop)Transpose to store (WBetaPROP)^t
  //NOTE the current will always have the pointers of the current WBeta where
  //beta is the current BETA and is an external object (reference)
  //(IMPORTANT this link the betas to the DPs and treatments)
  std::vector<std::unique_ptr<arma::Mat<double> >> &   RefVecPtrsAllTreatWBETAtransCurrent;
  std::vector<std::unique_ptr<arma::Mat<double> >>  VecPtrsAllTreatWBETAtransProposed;



  //Class for MH steps (Rylessone classes since we propose the nProgn betas simultaneously but
  //saparaterly for each category of the response )
  std::vector<MHProposal>  AllCatMHbetas;

  //Acceptance rates: is a vec of Ry-1 (since the BETA are updated joinly for each r)
  arma::Col<double> Acc;
  arma::Col<double> AccTotalprop;

  //Pointers and references to compute WtildeBetacurrent
  const arma::Mat<double> & RefWtildetrans;
  arma::Mat<double> & RefWBetaCurrtildetrans; //note that here need Wtildebeta need to be modifiable



public:
  //constructor Except
  GlobalParameters(Rcpp::List listpriors,
                   arma::uword const & ryless1,
                   arma::uword const & ntreat,
                   arma::uword const & nprognosticvar,
                   std::vector<std::unique_ptr<const arma::Mat<double> >>  const & PrognAllTreatPtr,
                   std::vector<std::unique_ptr<State>>  const & StateAllParametersExceptBeta,
                   std::vector<std::unique_ptr<arma::Mat<double> >>   & ExternalWBETACurrentTransAllTreatPtr ,
                   const arma::Mat<double> & refPostPredProgn,    arma::Mat<double> & refLinPredwbPostPred,
                   arma::Mat<double> Betainit, double lamdBeta):
  Rylessone(ryless1), nProgn(nprognosticvar),
  ntreatments(ntreat), //StateMVBPRtreat1(state),
  //MHbetas(arma::eye(nProgn,nProgn)*lamdBeta, 7 ),
  STATEalltreatmets_ptr(StateAllParametersExceptBeta),
  WtransAllTreatments_ptr(PrognAllTreatPtr),
  RefVecPtrsAllTreatWBETAtransCurrent(ExternalWBETACurrentTransAllTreatPtr),
  RefWBetaCurrtildetrans(refLinPredwbPostPred), RefWtildetrans(refPostPredProgn){
    //check if n treatments >0
    if(ntreatments<=0){
      throw std::runtime_error("n treatments must be > 0");

    }


    //check if the # prognostic variables is equal in all treatment: by default
    //nProgn is nrow of Wtrans^{treat=0}
    //also we save  the number of units in each treatment implied by prognostic variables
    arma::Col<arma::uword> NUnitsAllTreat(ntreatments);

    for (arma::uword nt = 0; nt < ntreatments; nt++) {
      NUnitsAllTreat(nt)=(*PrognAllTreatPtr.at(nt)).n_cols; //n units in each treatment for debug

      if( (*PrognAllTreatPtr.at(nt)).n_rows != nProgn){
        std::cout <<  "treatment nt " << nt << endl;
        throw std::runtime_error("N rows in Wtrans treatment nt != Nprong");

      }
    }

    //check if the predictive dataset the # prognostic variables is equal to
    if(RefWtildetrans.n_rows != nProgn){
      throw std::runtime_error("N rows in Wtildetrans != Nprong");

    }

    //check if StateAllParametersExceptBeta.nelem== Ntreat
    if(StateAllParametersExceptBeta.size() != ntreatments){
      throw std::runtime_error("Number of MVBPR not equal to Ntreat");
    }


    //read hyperparameters
    if (listpriors.containsElementNamed("HParbetas")) {
      // should be list of list
      Rcpp::List listHYbetas=listpriors["HParbetas"];
      //theta0
      if(listHYbetas.containsElementNamed("mean")){
        beta0=as<arma::Mat<double>>(listHYbetas["mean"]);
      }else{
        throw std::runtime_error("hyperparameter mean for beta's not given");
      }
      //scale
      if(listHYbetas.containsElementNamed("scale")){
        scale=as<arma::Mat<double>>(listHYbetas["scale"]);
      }else{
        throw std::runtime_error("hyperparameter scale for beta's not given");
      }
      //df
      if(listHYbetas.containsElementNamed("df")){
        df=as<arma::Mat<double>>(listHYbetas["df"]);
      }else{
        throw std::runtime_error("hyperparameter df for beta's not given");
      }

      //checks
      if(beta0.n_rows!=nProgn || beta0.n_cols!=Rylessone ){
        throw std::runtime_error("hyperparameter for betas: beta0.n_row!=nProgn || beta0.n_row!=ryless1");
      }
      if(scale.n_rows!=nProgn || scale.n_cols!=Rylessone ){
        throw std::runtime_error("hyperparameter for betas: scale.n_row!=nProgn || scale.n_row!=ryless1");
      }
      if(df.n_rows!=nProgn || df.n_cols!=Rylessone ){
        throw std::runtime_error("hyperparameter for betas: df.n_row!=nProgn || df.n_row!=ryless1");
      }

      //  if(arma::any(df < 2) || arma::any(scale < 0)){
      //  throw std::runtime_error("hyperparameter for betas: df |or| scale <0");
      // }
    }else{
      throw std::runtime_error("hyperparameters of betas not given");
    }

    //compute constant quantities for the prior
    //pre-compute normalizing constant and inv_dfscale
    inv_dfscale2=1.0 / (df % arma::square(scale));
    NegdfPlusOnedivTwo= -0.5*(df+1);
    //NorConst= arma::accu(arma::lgamma(0.5*(df+1))-arma::lgamma(0.5*df)-0.5*inv_dfscale2);
    logNorConst=arma::accu( arma::lgamma(0.5*(df+1))-arma::lgamma(0.5*df)-
      0.5*arma::log(arma::datum::pi*df) - arma::log(scale));




    //Initialize the Acceptance count
    Acc.resize(Rylessone);
    Acc.zeros();
    AccTotalprop.resize(Rylessone);
    AccTotalprop.zeros();


    //initialize beta
    BetasMatCurr=Betainit;
    BetasMatProp=Betainit;

    //Initialize MH classes
    AllCatMHbetas.reserve(Rylessone);

    for (arma::uword r = 0; r < Rylessone; r++) {
      AllCatMHbetas.emplace_back(arma::eye(nProgn,nProgn)*lamdBeta*1.0);
    }





    //initialize (Wbetaprop)^t  for all treatments
    VecPtrsAllTreatWBETAtransProposed.resize(ntreatments);

    //initialize as initial values we consider the Betainit
    for (arma::uword nt = 0; nt < ntreatments; nt++) {
      // NOTE: (W*BetaProp)^t =Beta_prop^t*Wtrans, and Wtrans is stored in the matrices accessible via pointers
      //   VecPtrsAllTreatWBETAtransCurrent.at(nt) = std::make_unique<arma::Mat<double>>(
      //            BetasMatCurr.t() * (*WtransAllTreatments_ptr.at(nt)) );

      VecPtrsAllTreatWBETAtransProposed.at(nt) = std::make_unique<arma::Mat<double>>(
        BetasMatProp.t() * (*WtransAllTreatments_ptr.at(nt)) );

      (*VecPtrsAllTreatWBETAtransProposed.at(nt)).fill(-5.0); //Here -5 anly for debug and this vector will be change by calling UPDATE





      //check if all the rows of Wbetas^{treat} are equal to Ry-1 (poposed)
      if((*VecPtrsAllTreatWBETAtransProposed.at(nt)).n_rows !=  Rylessone){
        std::cout <<  "treatment nt " << nt << endl;
        throw std::runtime_error(" nrows (WBeta)^trans for treatment nt != Ry-1");
      }

      //check if all the cols of Wbetas^{treat} are equal to NUnitsAllTreat (implies by W)
      //this is just for a pedantic check
      if((*VecPtrsAllTreatWBETAtransProposed.at(nt)).n_cols !=  NUnitsAllTreat.at(nt)){
        std::cout <<  "treatment nt " << nt << endl;

        throw std::runtime_error(" nrows (WBeta)^trans for treatment nt != N for treatment nt");
      }

      //check if the WBprop and WBcurr are same dimentions
      if( (*VecPtrsAllTreatWBETAtransProposed.at(nt)).n_cols !=  (*RefVecPtrsAllTreatWBETAtransCurrent.at(nt)).n_cols){
        std::cout <<  "treatment nt " << nt << endl;

        throw std::runtime_error(" ncol (WBetaProp)^trans for treatment nt != ncol (WBetaCurr)^trans for treatment nt");
      }
      if( (*VecPtrsAllTreatWBETAtransProposed.at(nt)).n_rows !=  (*RefVecPtrsAllTreatWBETAtransCurrent.at(nt)).n_rows){
        std::cout <<  "treatment nt " << nt << endl;

        throw std::runtime_error(" n_rows (WBetaProp)^trans for treatment nt != n_rows (WBetaCurr)^trans for treatment nt");
      }

    }


    cout << "-> Initialization  GlobalParameters class (class of the beta).... COMPLETE" <<endl;



  }
  //reset acc counter
  void ResetAccCounter(){
    Acc.zeros();
    AccTotalprop.zeros();
  }



  //add current beta to the accumulator to estimate the sample covariance
  void AccumulateBetasToEstimateRunningCov(){
    //initialize CubeAllcarAccumBetarBetartrans and Ma
    if(SampleSizeVariance==0){
      CubeAllcarAccumBetarBetartrans.emplace(nProgn,nProgn,Rylessone);
      MatrixAllCatAccumBetar.emplace(nProgn,Rylessone);
    }


    //Accumulate beta_r*beta_r^T (accumulate a matrix for each r)
    for( arma::uword r=0; r<Rylessone;r++){
      (*CubeAllcarAccumBetarBetartrans).slice(r)+=BetasMatCurr.col(r)*BetasMatCurr.col(r).t();

    }

    //Accumulate beta_r (accumulate a vector for each r, but we can directly sum Betacurr)
    (*MatrixAllCatAccumBetar)+=BetasMatCurr;

    SampleSizeVariance++;

  }


private:
  arma::mat MakeSPD(const arma::mat& A, double epsilon = 1e-8, bool add_jitter_if_needed = true) {
    // Step 1: Force symmetry
    arma::mat symA = 0.5 * (A + A.t());

    // Step 2: Eigendecomposition
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, symA);

    // Step 3: Clamp small/negative eigenvalues
    arma::vec clamped_eigval = arma::clamp(eigval, epsilon, arma::datum::inf);
    arma::mat repaired = eigvec * arma::diagmat(clamped_eigval) * eigvec.t();

    // Step 4: Optional jitter if smallest eigenvalue was below threshold
    if (add_jitter_if_needed && arma::min(eigval) < epsilon) {
      repaired += epsilon * arma::eye(A.n_rows, A.n_cols);
    }

    return repaired;
  }



  //compute Cov matrix using the accumulated betas
  arma::Cube<double> ComputeSampleCovMatrix(){
    arma::Cube<double> CovMathat(nProgn,nProgn,Rylessone);


    for( arma::uword r=0; r<Rylessone;r++){
      //Var_hat(X)=1/(P-1)  [Sum_i X*X^t - 1/n (Sum_i X_i)*(Sum_i X_i)^t]
      CovMathat.slice(r)=(*CubeAllcarAccumBetarBetartrans).slice(r);
      CovMathat.slice(r)-= (1/SampleSizeVariance)*(((*MatrixAllCatAccumBetar).col(r))*((*MatrixAllCatAccumBetar).col(r)).t());
      CovMathat.slice(r)=CovMathat.slice(r)/(SampleSizeVariance-1);

      //force symmetry
      CovMathat.slice(r) =MakeSPD(CovMathat.slice(r));
    }




    //reset Cube and matrix accumulator? mettere
    cout <<CovMathat<<endl;


    return CovMathat;
  }







  //return prior t-student computed in betas
  // as prior we assume independent beta_p,r~lst(beta0,s0r,df)
  double logPriorindtstud(const arma::Mat<double>  & betamat){

    arma::Mat<double> MatLogKern=NegdfPlusOnedivTwo % arma::log(1+ (arma::square(betamat-beta0) % inv_dfscale2));

    return (arma::accu(MatLogKern)+logNorConst);
  }

public:  
  //change the proposal variance using the sample covariance
  void ChangePropVariance(){
    
    auto EstimCov=ComputeSampleCovMatrix();
    
    for( arma::uword r=0; r<Rylessone;r++){
      AllCatMHbetas.at(r).ChangeScaleProposalUsingCovMat(EstimCov.slice(r));
      
    }
    
    
  }
  

  //Compute linear predictor WtildeBeta transposed
  void RecomputeWtildeBetaCurrPostPredtrans(){
    //This change the Matrix WtildeBetacurr [NPred x Ry-1], note that this matrix is
    //stored in In the main function not in the class, so can be linked to the STATE classes too

    RefWBetaCurrtildetrans=BetasMatCurr.t()*RefWtildetrans ;


  }

  //get acc
  arma::Mat<double> getAcc(){
    arma::Mat<double> MatAccTotProp(Acc.n_rows, 2);
    MatAccTotProp.col(0)=Acc;
    MatAccTotProp.col(1)=AccTotalprop;
    return MatAccTotProp;
  }

  //Update Beta
  void UPDATEBETAS(){

    //local variables to store the proposed logprior and loglikY associated to the BETATprop
    double logr, loglikYallTreatcurrent;//logpriorProp, loglikYallTreatProp, ;






    //compute current logprior   (we can avoid maybe...
    double logpriorcurrent=logPriorindtstud(BetasMatCurr);

    loglikYallTreatcurrent=0;
    //recover current loglikY all treatments (sum_trat loglikY_treat)
    for(auto& StateTthTreat : STATEalltreatmets_ptr){
      loglikYallTreatcurrent+=StateTthTreat->GetLogLikYallClusters();


    }

    //NOTA SI POTREBBE SOLAMENTE CONSIDERARE WBeta_r invece che BETAmat
    //ma per fare questo si deve cambiare logpriorprop e
    //        *VecPtrsAllTreatWBETAtransProposed.at(tr)=BetasMatProp.t() * (*WtransAllTreatments_ptr.at(tr));



    // for each beta_r [nProgn X 1]
    for(arma::uword r=0; r<Rylessone; r++){
      //  cout <<"  \n   ->updating cat "<<r  <<endl;

      //change r-th col of BetasMatProposed (that must, at this stage, to BetaMate current)
      //Propose new col beta_r, note betaMat_prop/curr [size: nProgn x Ry-1]
      BetasMatProp.col(r)=AllCatMHbetas.at(r).DrawMVT(BetasMatCurr.col(r));



      //compute propose logprior
      double logpriorProp=logPriorindtstud(BetasMatProp);

      //compute likY using the ``proposed'' (WBetaprop)^t for each treatment
      double loglikYallTreatProp=0;
      for(arma::uword tr=0; tr<ntreatments;tr++){


        //compute the new WBetaProp (trans) for t-th treatment using the new BetaProp
        *VecPtrsAllTreatWBETAtransProposed.at(tr)=BetasMatProp.t() * (*WtransAllTreatments_ptr.at(tr));

        //compute loglikY treatment t
        loglikYallTreatProp+=(*STATEalltreatmets_ptr.at(tr)).ComputeLogLikYallClustersGivenBetanew( *VecPtrsAllTreatWBETAtransProposed.at(tr));
      }


      //compute log ratio
      logr=logpriorProp+loglikYallTreatProp-logpriorcurrent-loglikYallTreatcurrent+
        AllCatMHbetas.at(r).logPropCorrFactor(BetasMatProp.col(r),BetasMatCurr.col(r));



      if(log(arma::randu())<logr){//accept proposed with prob min(1, exp(logr))
        // cout <<"-ACC-" <<endl;

        //if betaprop new accepted
        Acc(r)++;

        //change Betacurrent with the proposed, ie change the r-th col
        BetasMatCurr.col(r)=BetasMatProp.col(r);  //better swap matrix ??? more efficient

        //change current logprior and loglikY
        logpriorcurrent=logpriorProp;
        loglikYallTreatcurrent=loglikYallTreatProp;
        //      std::cout <<  "accetted ::rlog "<<logr  << "betmat"<<BetasMatCurr << endl;

        //Swap the pointers of (BETAcurrtrans*Wtrans)_treat t with the pointers of
        //(BETAProptrans*Wtrans)_treat t.
        //Afrer swap: the VecPtrsAllTreatWBETAtransCurrent store the pointers of the matrices previously
        //ponted by the pointes in   VecPtrsAllTreatWBETAtransProposed
        //Instead the VecPtrsAllTreatWBETAtransProposed store the matrices previously pointed by the
        //pointers inVecPtrsAllTreatWBETAtransCurrent ,
        // but is not a problem since at the next step will be recomputed
        std::swap(RefVecPtrsAllTreatWBETAtransCurrent, VecPtrsAllTreatWBETAtransProposed);





      }else{
        //if betaprop not accepted
        //reset the BetasMatProposed==BetasMat since at the next iteration
        //the r-th+1 col will be change, but the others need to be equal to beta current
        BetasMatProp.col(r)=BetasMatCurr.col(r);

        //note in this case loglikYallTreatcurrent and log logpriorcurrent ar OK
        //while the logpriorProp loglikYallTreatProp are no longer valid since
        //after   BetasMatProp.col(r)=BetasMatCurr.col(r), the two martix are equal
        //this is not a proble since at the nex step the proposed quantities will be
        //change and loglik and log prop follows in turn

        //    std::cout <<  "NOT accepted ::rlog "<<logr  << "betmat"<<BetasMatCurr << endl;

      }

    }
    /*
     cout <<"\n ---end of BETA updating:---" <<endl;
     cout <<"pointer vec WB current "<< RefVecPtrsAllTreatWBETAtransCurrent.at(0).get() <<endl;
     cout <<"pointer vec WB proposed "<<VecPtrsAllTreatWBETAtransProposed.at(0).get()  <<endl;


     cout <<"current loglikY "<< loglikYallTreatcurrent <<endl;

     cout <<"current loglikY recomputed "<< loglikYallTreatcurrent <<endl;

     cout <<"\n \n "<<endl;
     */

  }

  //save betas
  void SaveGlobalParms(arma::cube & Bmatpost, arma::uword mcmcind){
    Bmatpost.slice(mcmcind)=BetasMatCurr;

  }





};






//Sample from the Y model: Y~cat(eta) with eta_r \propto exp(theta_kr+Wbeta_r)
arma::Col<arma::uword> SampleResponse(//arma::Col<double> const & Thetak,
    //        arma::uword const & npred,
    arma::Mat<double>  & UnnormLogProbs){
  //note that

  // Append zero log-prob column for last category
  //  UnnormLogProbs = arma::join_horiz(UnnormLogProbs, arma::zeros(npred));

  // Generate Gumbel noise: G = -log(-log(U)), U ~ Uniform(0,1) and
  //add to the matrix
  UnnormLogProbs+= -arma::log(-arma::log(arma::randu<arma::mat>(arma::size(UnnormLogProbs))));


  return arma::index_max(UnnormLogProbs, 1);

}

 



arma::Mat<arma::uword> InitZmatrix(arma::uword n, arma::uword group, arma::uword k = 2) {
  arma::Mat<arma::uword> out(n, k);

  for (arma::uword j = 0; j < k; ++j) {
    arma::Col<arma::uword> col(n);
    // First, assign each group label once (for coverage)
    for (arma::uword g = 0; g < group; ++g) {
      col(g) = g + 1;
    }
    // Assign the remaining (n - group) randomly
    col.subvec(group, n - 1) = arma::randi<arma::Col<arma::uword>>(n - group, arma::distr_param(1, group));
    // Shuffle the column so guaranteed labels are mixed
    col = arma::shuffle(col);

    out.col(j) = col;
  }

  return out;
}
 
 // [[Rcpp::export]]
 Rcpp::List MVBPRmultitreat(List DataList,
                            int M,
                            double lambdtheta, //  variance proposal theta
                            double lambdbeta, // variance proposal betas
                            int numbtreatments, // 
                            int numbprognosticvar,
                            NumericVector maximumcatDiscrete, //
                            int NumberofView, //
                            int TypeofXmodel, //OK MUST BE ALWAYS set to 1
                            int Maxcatresponse, //
                            List ListPrior, //OK  ,
                            int Ninitclusters,
                            int InitGammas, //if -2 matrix ,-1 random otherwise all D categorical vars allocated to 0 1 2 .. Nview-1
                            NumericMatrix InitGammasmat, //if InitGammas==-2 
                            NumericVector AlphaDPsinit,
                            NumericMatrix BetaMatinit,
                            int MCMCFinalSamplesize,
                            int Burnin=0,
                            int Thinning=1,
                            int ReturnBurnin=0){ //if 0 no, 1 yes all burn in, 2 only after change proposal variance  
   
 

   //check positiva value MCMC parameters
   if(Thinning<1 || MCMCFinalSamplesize<0 || Burnin<0 ){
     throw std::runtime_error("MVBPRmultitreat: Thinning<1 || MCMCFinalSamplesize<0 || Burnin<0");
     
   }
   if(numbprognosticvar<=0 ){
     throw std::runtime_error("MVBPRmultitreat: numbprognosticvar must be > 0");
     
   }
   
   
   if(InitGammas< -2 || InitGammas>=NumberofView){
     throw std::runtime_error("MVBPRmultitreat: InitGammas< -1 || InitGammas>=NumberofView");
     
   }
   
  //create a Rcpp list to save the names of Xdisc Xcont and W
  Rcpp::List ListVariablesNames;
   
   //recompute MCMCFinalSamplesize
   int MCMCsaveburnin;
   if(ReturnBurnin==1){
     
     MCMCsaveburnin=Burnin*3;
     
   } else if(ReturnBurnin==2){
     
     MCMCsaveburnin=Burnin;
     
   }else{
     
     MCMCsaveburnin=0;
     
   }
   MCMCFinalSamplesize+=MCMCsaveburnin;
   
   
   //-----------------global non-treatment data--------------------------------
   //N categorical predictive variables
   const  arma::uword Dvar=static_cast<arma::uword>(maximumcatDiscrete.length());
   
   //ATTENZIONE conversion from ‘R_xlen_t’ {aka ‘long int’} to ‘arma::uword’ {aka ‘unsigned int’} may change value [-Wconversion]
   
   //max categories X_d discrete
   const arma::Col<arma::uword> armamaxd = as<arma::Col<arma::uword>>(maximumcatDiscrete);
   //categories response variable less 1
   const arma::uword armamaxless1Y=static_cast<arma::uword>(Maxcatresponse-1);
   
   //convert to arma
   arma::uword NprognVar=static_cast<arma::uword>(numbprognosticvar);
   arma::uword Ntreatments=static_cast<arma::uword>(numbtreatments);
   
   arma::uword NVinit=static_cast<arma::uword>(NumberofView);
   //-----------------Treatment specific data--------------------------------
   
   //Convert data in ARMA structures  (DA MODIFICARE)
   //const arma::Mat<arma::uword> armaDatad = as<arma::Mat<arma::uword>>(Datad);
   //const arma::Mat<double> armaDatac = as<arma::Mat<double>>(Datac);
   //const arma::Col<arma::uword> armaY = as<arma::Col<arma::uword>>(Responsevec);
   
   //prognostic variables one matrix for each treatment (each matrix has dim [Nprogn X Nstate])
   /* std::vector<std::shared_ptr<const arma::Mat<double> >> PrognAllTreatmetsTrans;
    PrognAllTreatmetsTrans.resize(Ntreatments);
    for (arma::uword nt = 0; nt < Ntreatments; nt++) {
    PrognAllTreatmetsTrans.at(nt)=std::make_shared<arma::Mat<double>>(as<arma::Mat<double>>(Progn).t());
    }*/
   
   //-----------------Treatment specific data--------------------------------
   //save data (cont-disc-progn) in a const vector of unique pointers. each pointer
   //point to a const  const arma::Mat<double>
   
   //init temp vector of pointers mat double (used for cont)
   std::vector<std::unique_ptr<const arma::Mat<double>>> tempPtrVectCONT;
   tempPtrVectCONT.resize(Ntreatments);
   //init temp vector of pointers mat double (used for progn)
   std::vector<std::unique_ptr<const arma::Mat<double>>> tempPtrVectPROGN;
   tempPtrVectPROGN.resize(Ntreatments);
   //init temp vector of pointers arma::Mat<arma::uword> (used for disc)
   std::vector<std::unique_ptr<const arma::Mat<arma::uword>>> tempPtrVectDISC;
   tempPtrVectDISC.resize(Ntreatments);
   //init temp vector of pointers arma::Col<arma::uword> (used for response)
   std::vector<std::unique_ptr<const arma::Col<arma::uword>>> tempPtrVectRESP;
   tempPtrVectRESP.resize(Ntreatments);
   
   
   
   //populate the vectors by creating (converting) the ARMA data matrix for each treatments
   //each component (X_d^t, W^t,y^t) is convetred from Rcpp Matrix/vector (stored in the list) to an arma equivalent
   for (arma::uword nt = 0; nt < Ntreatments; nt++) {//MODIFICARE PER LEGGERE TREARMENT SPECIFIC DATA
     //read data to fit the model from DataList
     Rcpp::List fitdata = DataList["fitdata"]; //this line could be place outside the for but in this way we free memory after reading is completed
     Rcpp::List Treatment_tDataFit = fitdata[nt];                         // treatment specific list: X Y W
     
     //discrete variables
     Rcpp::NumericMatrix Datad = Treatment_tDataFit["XDiscPred"];  //Rcpp Matrix X discrete
     
     //Rcpp::CharacterVector colnames = Rcpp::clone(Datad.attr("dimnames").get());
     
     
     tempPtrVectDISC.at(nt)=std::unique_ptr<const arma::Mat<arma::uword>>(
       new  arma::Mat<arma::uword>(as<arma::Mat<arma::uword>>(Datad)));
     
     //continuous variables NOT IMPLEMENTED
     Rcpp::NumericMatrix Datac = Treatment_tDataFit["XContPred"];  //Rcpp Matrix X cont. Placeholder [n x 0] matrix in this implementation
     
     //cout << Rcpp::List(Datad.attr("dimnames"))[1]. <<endl;
   //  cout << Datad.ncol()  <<Datad.nrow()  <<endl;
     
     tempPtrVectCONT.at(nt)=std::unique_ptr<const arma::Mat<double>>(
       new  arma::Mat<double>(as<arma::Mat<double>>(Datac)));
     
     //prognostic variables (FOR THE PROGNOSTIC VARIBLES we consider the transpose... each matrix has dim [Nprogn X Nstate])
     Rcpp::NumericMatrix Progn = Treatment_tDataFit["WProgn"];  //Rcpp Matrix W
     
     tempPtrVectPROGN.at(nt)=std::unique_ptr<const arma::Mat<double>>(
       new arma::Mat<double>(as<arma::Mat<double>>(Progn).t()*static_cast<double>(1.0)));
     
     //save variables name
     Rcpp::CharacterVector colnamesxd = 
           Datad.ncol() > 0 
         ? Rcpp::as<Rcpp::CharacterVector>(Rcpp::List(Datad.attr("dimnames"))[1])
           : Rcpp::CharacterVector();  
         
     
     Rcpp::CharacterVector colnamesxc = 
           Datac.ncol() > 0  
         ? Rcpp::as<Rcpp::CharacterVector>(Rcpp::List(Datac.attr("dimnames"))[1])
           : Rcpp::CharacterVector();
     
     Rcpp::CharacterVector colnamesw =   
           Progn.ncol() > 0  
          ? Rcpp::as<Rcpp::CharacterVector>(Rcpp::List(Progn.attr("dimnames"))[1])
            : Rcpp::CharacterVector();  
     
     ListVariablesNames["XdiscNames"]=Rcpp::clone(colnamesxd);
     ListVariablesNames["XcontNames"]=Rcpp::clone(colnamesxc);
     ListVariablesNames["XprogNames"]=Rcpp::clone(colnamesw);
     
     //Response
     Rcpp::NumericVector Responsevec=Treatment_tDataFit["Y"];
     
     tempPtrVectRESP.at(nt)=std::unique_ptr<const arma::Col<arma::uword>>(
       new  arma::Col<arma::uword>(as<arma::Col<arma::uword>>(Responsevec)));
     
     
     
     
     //Check dim dataset not that further check will be performed by the class initializator  
     if((*tempPtrVectDISC.at(nt)).n_rows != (*tempPtrVectPROGN.at(nt)).n_cols){
       std::cout << "nt="<< nt << endl;
       throw std::runtime_error("MVBPRmultitreat: in trt nt, nrow of Disc Predictive X and Prognostic W are NOT equal");
       
     }
     if( (*tempPtrVectDISC.at(nt)).n_cols != Dvar){  
       std::cout << "nt="<< nt << endl;
       throw std::runtime_error("MVBPRmultitreat: ncol of Disc Predictive X not equal to number of predictive X in the model, ie maximumcatDiscrete.length");
       
     }
     
     if( (*tempPtrVectPROGN.at(nt)).n_rows != NprognVar){
       std::cout << "nt="<< nt << endl;
       throw std::runtime_error("MVBPRmultitreat: N of Prognostic covariate not equal to number of prognostic variables in the model, ie numbprognosticvar");
       
     }
     
   }
   
   //MOVE THESE POINTERS IN A CONST VECTORS
   const std::vector<std::unique_ptr<const arma::Mat<arma::uword>>> AllTreatDataDiscVectOfPtr(std::move(tempPtrVectDISC));
   const std::vector<std::unique_ptr<const arma::Mat<double>>> AllTreatDataContVectOfPtr(std::move(tempPtrVectCONT));
   const std::vector<std::unique_ptr<const arma::Mat<double>>> AllTreatDataProgTransVectOfPtr(std::move(tempPtrVectPROGN));
   const std::vector<std::unique_ptr<const arma::Col<arma::uword>>> AllTreatRespVectOfPtr(std::move(tempPtrVectRESP));
   
   
   
   //-----------------Global parameters initialization------------------------
   Rcpp::List PostPreddata = DataList["postpredictivedata"];
   //Rcpp::NumericMatrix Datacp = PostPreddata["XDiscPred"];
   //Rcpp::NumericMatrix Prognp = PostPreddata["WProgn"];  //Rcpp Matrix W
   
   //pointers Matrix X and W for of the predictive units [NpostPred x Dvar]
   const std::unique_ptr<const arma::Mat<arma::uword>> PtrPostPredDatac(new arma::Mat<arma::uword>(as<arma::Mat<arma::uword>>(PostPreddata["XDiscPred"])));
   
   //Pointer prognostic variables  transposed, matrix [Nprogn x NpostPred]
   const std::unique_ptr<const arma::Mat<double>> PtrPostPredProgntrans(new arma::Mat<double>(as<arma::Mat<double>>(PostPreddata["WProgn"]).t()));
   
   
   
   
   
   //Check on post predictive dataset
   if((*PtrPostPredDatac).n_rows != (*PtrPostPredProgntrans).n_cols){
     throw std::runtime_error("error: n observations of the Post Predictive X covariate not equal to the Post Predictive Prognostic");
     
   }
   if((*PtrPostPredDatac).n_cols != Dvar){ 
     cout << "ncol x tilde" << (*PtrPostPredDatac).n_cols << "while Dvar " <<Dvar <<endl;
     throw std::runtime_error("error: N colss Post Predictive X covariate not equal to number of predictive X in the model");
     
   }
   
   if((*PtrPostPredProgntrans).n_rows != NprognVar){
     
     
     cout << NprognVar <<"-" <<(*PtrPostPredProgntrans).n_rows <<endl;
     throw std::runtime_error("error: N rows Post Predictive Prognostic covariate not equal to number of prognostic variables in the model");
     
   }
   
   
   
   //Reference (current) Linear predictor WtildeBetacurrent Posterior predictive [Ry-1 x NpostPred]
   arma::Mat<double> WPostPredBetacurrtrans = arma::Mat<double>( armamaxless1Y, (*PtrPostPredDatac).n_rows);
   WPostPredBetacurrtrans.zeros();
   
   
   
   // THIS LINK THE MVBPR CLASSES WITH THE GLOBAL PARAMETERS CLASS
   
   //initialize (Wbetaprop)^t  for all treatments: TWO cases: 1 no BETA or 2 With beta
   std::vector<std::unique_ptr<arma::Mat<double> >>  VecPtrsAllTreatWBETAtransCurrent;
   
   VecPtrsAllTreatWBETAtransCurrent.resize(Ntreatments);
   
   
   
   //we consider the 1st treatment, all treatment must have the same nprogn
   // const arma::uword Nprognosticvariables=(*AllTreatDataProgTransVectOfPtr.at(0)).n_cols;
   
   //cout << "Nprognosticvariables" <<Nprognosticvariables<<endl;
   //define matrix beta (used only if Npron!=0)
   arma::Mat<double> Betainitmatarma;
   
   //initialize as initial values we consider the Betainit
   for (arma::uword nt = 0; nt < Ntreatments; nt++) {
     //IF NprognVar!=0
     if(NprognVar!=0){
       Betainitmatarma= as<arma::Mat<double>>(BetaMatinit);
       //check if the matrix BETAinit dims is coherent with Ry-1, and n prognostic var
       if(Betainitmatarma.n_cols!=armamaxless1Y || Betainitmatarma.n_rows!= NprognVar){
         throw std::runtime_error("error: Betainitmatarma.n_cols!=NprognVar ||Betainitmatarma.n_rows!= armamaxless1Y");
         
       }
       //check if the matrix (*AllTreatDataProgTransVectOfPtr.at(nt)) dims is coherent with Ry-1, and n prognostic var
       if( (*AllTreatDataProgTransVectOfPtr.at(nt)).n_rows!= NprognVar){
         throw std::runtime_error("error: (*AllTreatDataProgTransVectOfPtr.at(nt)).n_rows!= NprognVar");
         
       }
       // NOTE: (W*BetaProp)^t =Beta_prop^t*Wtrans, and Wtrans is stored in the matrices accessible via pointers
       //   VecPtrsAllTreatWBETAtransCurrent.at(nt) = std::make_unique<arma::Mat<double>>(
       //            BetasMatCurr.t() * (*WtransAllTreatments_ptr.at(nt)) );
       
       VecPtrsAllTreatWBETAtransCurrent.at(nt) = std::make_unique<arma::Mat<double>>(
         Betainitmatarma.t() * (*AllTreatDataProgTransVectOfPtr.at(nt)) );
       
     }else{
       //IF there are no prognostic variables
       //VecPtrsAllTreatWBETAtransCurrent.at(nt) = std::make_unique<arma::Mat<double>>(
       //     arma::Mat<double>( arma::Mat<double>(armamaxless1Y,(*AllTreatDataProgTransVectOfPtr.at(nt)).n_cols )));
       
       
       
       //for debug
       Betainitmatarma= as<arma::Mat<double>>(BetaMatinit);
       
       VecPtrsAllTreatWBETAtransCurrent.at(nt) = std::make_unique<arma::Mat<double>>(
         Betainitmatarma.t() * (*AllTreatDataProgTransVectOfPtr.at(nt)) );
       
     }
     
     cout << "Treat" << nt << "dim " << (*VecPtrsAllTreatWBETAtransCurrent.at(nt)).n_rows << "," <<
       (*VecPtrsAllTreatWBETAtransCurrent.at(nt)).n_cols <<endl;
     
   }
   
   //Initialize list to save inital gammas and z
   Rcpp::List InitialVal;
   
   
   //-----------------Treatment specific initialization------------------------
   //  (DA MODIFICARE: devono essere treatment specific Z and alphas)
   //  arma::Col<arma::uword> Garma = as<arma::Col<arma::uword>>(Gammasinit);
   // arma::Mat<arma::uword> Zarma = as<arma::Mat<arma::uword>>(Zmatinit); //MODIFICARE DEVE ESSERE TREATMENT SPECIFIC
   arma::Col<double>  alphaDPsarma =as<arma::Col<double>>(AlphaDPsinit);
   
   
   
   arma::uword M8nealarma=static_cast<arma::uword>(M);
   
   
   //-----------------CLASS initializations ------------------------
   std::vector<std::unique_ptr<State>> tempMVBPRtempvec;  //temp  Vector to store unique pointers STATE MVBPR
   tempMVBPRtempvec.resize(Ntreatments);
   
   
   
   
   cout << " ____________________________________________________________________"<<endl;
   cout << "|       Multi-Treatment Multi-view Bayesian Profile Regression       |"<<endl;
   cout << "|____________________________________________________________________|"<<endl;
   cout << "   Number of treatments: " << Ntreatments  <<endl;
   cout << "   Number of Discrete predictive variables: " << Dvar  <<endl;
   cout << "   Number of prognostic variables: " << NprognVar  <<endl;
   
   
   
   //Initialize Class MVBPR for each treatments: saved in STATEalltreatmets as pointer
   for (arma::uword nt = 0; nt < Ntreatments; nt++) {
     cout<< " "<<endl;
     cout << "-------------------------"   <<endl;
     cout << "Treatment: " << nt  <<endl;
     cout << "-------------------------"   <<endl;
     cout << "  Number of observations:" << (*AllTreatDataDiscVectOfPtr.at(nt)).n_rows   <<endl;
     cout << "  Number of Views:" << NVinit    <<endl;
     
     //init cluster allocations
     auto Zarma=InitZmatrix((*AllTreatDataDiscVectOfPtr.at(nt)).n_rows, static_cast<arma::uword>(Ninitclusters)-1, NVinit-1); //arma::uword n, arma::uword group, arma::uword k = 2
     
     //init gammas
     arma::Col<arma::uword> Garma;
     if(InitGammas == -1){//random
       Garma= arma::randi<arma::Col<arma::uword>>(Dvar, arma::distr_param(0, NVinit-1));
     }else{
       if(InitGammas == -2){
         auto Gmat=as<arma::Mat<arma::uword>>(InitGammasmat);
         
         Garma=Gmat.col(nt);
         
       }else{
         Garma=arma::ones<arma::Col<arma::uword>>(Dvar) * InitGammas;
         
       }
     }
     
     
     std::string InitialValz = "trtt" + std::to_string(nt + 1)+ "_Z";
     std::string InitialValg = "trtt" + std::to_string(nt + 1)+ "_Gamma";
     
     InitialVal[InitialValz]=Zarma;
     InitialVal[InitialValg]=Garma;
     //  cout << Garma <<endl;
     //Garma(90)=1;
     //Garma(91)=1;Garma(92)=1;Garma(93)=1;Garma(94)=1;Garma(95)=1;Garma(96)=1;Garma(97)=1;
     tempMVBPRtempvec.at(nt)=std::make_unique<State>(Zarma, Garma, alphaDPsarma, //initial ZZ,Gammas,alpha_v
                         M8nealarma, //parameter m alghorithm 8 neal
                         (*AllTreatDataDiscVectOfPtr.at(nt)),//Data treatment specific DISCRETE
                         (*AllTreatDataContVectOfPtr.at(nt)),//Data treatment specific CONT (NOT IMPLEMENTED)
                         (*AllTreatRespVectOfPtr.at(nt)), //Data treatment specific RESPONSE
                         TypeofXmodel,              //global across treatment
                         NVinit,
                         armamaxd,
                         Dvar,
                         armamaxless1Y,
                         VecPtrsAllTreatWBETAtransCurrent,
                         nt,
                         ListPrior,
                         lambdtheta,
                         (*PtrPostPredDatac),
                         WPostPredBetacurrtrans);
   }
   
   //Convert tempMVBPRtempvec in constant vector so the pointer are immutable
   const std::vector<std::unique_ptr<State>> STATEMVBPRAllTreatVectOfPtr(std::move(tempMVBPRtempvec));
   
   std::optional<GlobalParameters> BETACLASS;
   //conditional initialization: in there are prognostc variables
   if(NprognVar!=0){
     cout << "-------------------------"   <<endl;
     cout << "Global Parameters class  "  <<endl;
     cout << "-------------------------"   <<endl;
     BETACLASS.emplace(ListPrior,
                       armamaxless1Y,
                       Ntreatments,
                       NprognVar,
                       AllTreatDataProgTransVectOfPtr,
                       STATEMVBPRAllTreatVectOfPtr,
                       VecPtrsAllTreatWBETAtransCurrent ,
                       *(PtrPostPredProgntrans), WPostPredBetacurrtrans,
                       Betainitmatarma, lambdbeta  );
   }
   
   
   
   
   
   
   //-----------------initializations structure to return------------------------
   
   //Initialize data structures to save the posterior drawn of Z_1, Z_2,...,Z_v, gammas for all treatments
   Rcpp::List LISTallTREATMENTS;
   //each element of the vector is a cube, each cube save the allocation variables
   //->ZZ_posteriorAllTreatments [N x nnView x MCMC]: row represent the i-th units, col represent the v-th view,
   // and the 3rd index represent the s-th saved draw from the posterior: Z_{v,i}^[s]
   //->Gammas_posteriorAllTreatments [D x MCMC]: row represent the i-th units, col represent
   // the s-th saved draw from the posterior
   //->Alpha_posterior [nnView x MCMC]: row represent the i-th units, col represent
   // the s-th saved draw from the posterior
   //->Theta_usAllTreatments [Ry-1 x N x MCMC]: row represent the r-th parameter (with r=1...Ry-1)
   //col represent the i-th unit, 3rd index the mcmc simulation
   std::vector<std::shared_ptr<arma::ucube>> ZZ_posteriorAllTreatments;
   std::vector<std::shared_ptr<arma::umat>> Gammas_posteriorAllTreatments;
   std::vector<std::shared_ptr<arma::mat>> Alpha_posteriorAllTreatments;//(alphaDPsarma.n_elem, MCMCsim)
   std::vector<std::shared_ptr<arma::cube>> Theta_usAllTreatments;
   
   
   
   std::vector<std::shared_ptr<arma::umat>> PostPredictive_posteriorAllTreatments;//(alphaDPsarma.n_elem, MCMCsim)
   
  // std::vector<std::shared_ptr<arma::cube>> PostPredictive_DEBUG;//(alphaDPsarma.n_elem, MCMCsim)
   
   
   //Acceptance counters
   
   // if(BETACLASS){
   
   std::shared_ptr<arma::cube> Betas_posterior=std::make_shared<arma::cube>(NprognVar, armamaxless1Y, MCMCFinalSamplesize);
   //}
   
   
   ZZ_posteriorAllTreatments.resize(Ntreatments);
   Gammas_posteriorAllTreatments.resize(Ntreatments);
   Alpha_posteriorAllTreatments.resize(Ntreatments);
   Theta_usAllTreatments.resize(Ntreatments);
   PostPredictive_posteriorAllTreatments.resize(Ntreatments);
   
  // PostPredictive_DEBUG.resize(Ntreatments);
   
   //DA MODIFICARE PER Zarma.n_rows TREATMENT SPECIFIC. FATTO???!
   for (arma::uword nt = 0; nt < Ntreatments; nt++) {
     
     ZZ_posteriorAllTreatments.at(nt)=std::make_shared<arma::ucube>((*AllTreatDataDiscVectOfPtr.at(nt)).n_rows,
                                  NVinit-1,
                                  MCMCFinalSamplesize);
     Gammas_posteriorAllTreatments.at(nt)=make_shared<arma::umat>(Dvar, MCMCFinalSamplesize);
     Alpha_posteriorAllTreatments.at(nt)=make_shared<arma::mat>(alphaDPsarma.n_elem, MCMCFinalSamplesize);
     Theta_usAllTreatments.at(nt)=make_shared<arma::cube>(armamaxless1Y, (*AllTreatDataDiscVectOfPtr.at(nt)).n_rows, MCMCFinalSamplesize);
     
     PostPredictive_posteriorAllTreatments.at(nt)=make_shared<arma::umat>(WPostPredBetacurrtrans.n_cols,
                                              MCMCFinalSamplesize);
     
     cout <<(*PostPredictive_posteriorAllTreatments.at(nt)).n_cols <<"-"<<(*PostPredictive_posteriorAllTreatments.at(nt)).n_rows <<endl;
     
     //    PostPredictive_DEBUG.at(nt)=make_shared<arma::cube>(WPostPredBetacurrtrans.n_cols, armamaxless1Y+1,
     //                          MCMCFinalSamplesize);
     
     
   }
   
   // arma::vec DEBUGlogliky(MCMCFinalSamplesize);
   
   
   //auxiliary variables for MCMC sampler
   arma::uword indextosave=0;
   int saveCounter=0;
   
 //iteration after burn in
   int MCMCTotIter= (MCMCFinalSamplesize-MCMCsaveburnin)*Thinning;//+Burnin) ;
   
   
   //auc variables to print   
   const uint64_t interval = (MCMCTotIter+Burnin*3) / 10;    
   uint64_t counterinterval = interval;                  
   uint64_t printtenpercent = 10;                        
   
   
   
   
   //#########################################
   //SAMPLER
   cout << "\n \n "     <<endl;
   cout << "_____________________________________ "     <<endl;
   cout << "____________START SAMPLER____________ "     <<endl;
   cout << "Burnin: "<< Burnin*3     <<endl;
   cout << "MCMC total iteration: "<<MCMCTotIter     <<endl;
   cout << "Thinning: " <<Thinning <<endl;
   cout << "MCMC Final Sample size: "<<MCMCFinalSamplesize     <<endl;
   cout << "\n  "     <<endl;
   
   
   
   
   
   std::clock_t start = std::clock();
   //----------------BURNIN---------------- not saving
   for(int s=0; s<Burnin*3 ;s++){
     //UPDATE MVBPR's all treatments (except Betas)
     for(auto& state : STATEMVBPRAllTreatVectOfPtr){
    
       state->SAMPLEROneRun();
       
     }
     
     //UPDATE BETAS
     if(NprognVar!=0){
       BETACLASS ->  UPDATEBETAS(); //UPDATEBETAS();
       
       
       //second part of burnin [Burnin, 2*Burnin] estimate var proposal
       if(s>=Burnin && s<Burnin*2){
         
         BETACLASS ->AccumulateBetasToEstimateRunningCov();
       }
       
      }
     
     
     

     if(ReturnBurnin==1 ||  (ReturnBurnin==2 && s>=Burnin*2 )  ){
       //Computing post predictive (updeting Wtildebeta with current) THIS STEP IS PARAMOUNT otherwise the Beta would not be considered
       BETACLASS->RecomputeWtildeBetaCurrPostPredtrans();
       
       //save thetak's gammas alphas etc+ post pred, loop over all treatments 
       for(arma::uword tr=0; tr<Ntreatments;tr++){
         
      /*   void SaveParamsandPostPred(arma::ucube & Zsave,
                                    arma::umat & GammasSave,
                                    arma::cube & Thetaus,
                                    arma::mat & Alphasave,
                                    arma::umat   & YPostpred,
                                    arma::uword mcmcind)*/
         
         STATEMVBPRAllTreatVectOfPtr.at(tr) ->SaveParamsandPostPred((*ZZ_posteriorAllTreatments.at(tr)),
                                                                    (*Gammas_posteriorAllTreatments.at(tr)),
                                                                    (*Theta_usAllTreatments.at(tr)),
                                                                    (*Alpha_posteriorAllTreatments.at(tr)),
                                                                    (*PostPredictive_posteriorAllTreatments.at(tr)),
                                                                    indextosave);
       }
       
       
       if(NprognVar!=0){
         BETACLASS ->SaveGlobalParms(*(Betas_posterior), indextosave);
         
       }
       
       //
       indextosave++;
     }
     
     //print percent completed
     if (--counterinterval == 0) {
       Rcpp::Rcout << "\r"
                   << printtenpercent
                   << "% complete";
       Rcpp::Rcout.flush();
       
       // std::cout << " \r" 
       //           << printtenpercent 
       //           << "% complete" 
       //           << std::flush;
       
       printtenpercent += 10;
       counterinterval = interval;
     }
   }
   
   
   
   
   
   //reset counters acceptance
   if(Burnin>0){
     for(auto& state : STATEMVBPRAllTreatVectOfPtr){
       state->ResetAccCounter();
     }
     if(NprognVar!=0){
       BETACLASS->ResetAccCounter();
       
       //set estimated var proposal
       BETACLASS ->ChangePropVariance();
     }
   }
   
   cout << " -Burnin Complete \n"     <<endl;
   
   
   
   //----------------MAIN (saving)---------------
   for(int s=0; s<MCMCTotIter;s++){//loop MCMC
     //UPDATE MVBPR's all treatments (except Betas)
     //loop over all treatments (update all parameters of MVBPR associated to treatment t: Zt,v_i's, Theta_k's=t)
     for(auto& state : STATEMVBPRAllTreatVectOfPtr){
       
       state->SAMPLEROneRun();
       
     }
     
     
     //UPDATE BETAS
     if(NprognVar!=0){
       BETACLASS->UPDATEBETAS();
     }
     
     
     
     //SAVE RESULTS (according to the Thinning) and Computing the Post predictive
     //increase saveCounter
     saveCounter++;
     
     if(saveCounter==Thinning   ){
       
       //Computing post predictive (updeting Wtildebeta with current) THIS STEP IS PARAMOUNT otherwise the Beta would not be considered
       BETACLASS->RecomputeWtildeBetaCurrPostPredtrans();
       
       //save thetak's gammas alphas etc+ post pred, loop over all treatments 
       for(arma::uword tr=0; tr<Ntreatments;tr++){
         
         /*   void SaveParamsandPostPred(arma::ucube & Zsave,
          arma::umat & GammasSave,
          arma::cube & Thetaus,
          arma::mat & Alphasave,
          arma::umat   & YPostpred,
          arma::uword mcmcind)*/
         
         STATEMVBPRAllTreatVectOfPtr.at(tr) ->SaveParamsandPostPred((*ZZ_posteriorAllTreatments.at(tr)),
                                        (*Gammas_posteriorAllTreatments.at(tr)),
                                        (*Theta_usAllTreatments.at(tr)),
                                        (*Alpha_posteriorAllTreatments.at(tr)),
                                        (*PostPredictive_posteriorAllTreatments.at(tr)),
                                        indextosave);
       }
       
       
       if(NprognVar!=0){
         BETACLASS ->SaveGlobalParms(*(Betas_posterior), indextosave);
         
       }
       
       
       //
       indextosave++;
       saveCounter=0;
       
     }
     
     
     //print percent completed
     if (--counterinterval == 0) {
       Rcpp::Rcout << "\r"
                   << printtenpercent
                   << "% complete";
       Rcpp::Rcout.flush();
       
      // std::cout << " \r" 
      //           << printtenpercent 
      //           << "% complete" 
      //           << std::flush;
       
       printtenpercent += 10;
       counterinterval = interval;
     }
     
   }//END SAMPLER
   cout << " \n "     <<endl;
   
   
   
   std::clock_t end = std::clock();
   double elapsed = static_cast<double>(end - start) / CLOCKS_PER_SEC;
   std::cout << "CPU Time: " << elapsed << " seconds" << std::endl;
   cout << "_____________________________________ "     <<endl;
   
   
   
   
   
   //
   //LIST to ruturn
   for(arma::uword tr=0; tr<Ntreatments;tr++){
     //convert Accthetas in vector to export in R
     
     const auto AccTheta_umap=(*STATEMVBPRAllTreatVectOfPtr.at(tr)).getAccThetak();//.AccThetak;
     
     Rcpp::NumericMatrix MatAccTheta_k(2, static_cast<int>(AccTheta_umap.size()));
     
     
     int j=0;
     for(const auto & accThetak: AccTheta_umap){
       //index of theta_k, ie k
       MatAccTheta_k(0, j)=accThetak.first;
       //accepted
       MatAccTheta_k(1,j)=accThetak.second(0);///accThetak.second(1);
       j++;
     }
     
     Rcpp::List ListToReturn;
     
     ListToReturn["ZZallviews"]=(*ZZ_posteriorAllTreatments.at(tr));
     ListToReturn["Gammas"]=(*Gammas_posteriorAllTreatments.at(tr));
     ListToReturn["alpha"]=(*Alpha_posteriorAllTreatments.at(tr));
     ListToReturn["Thetaus"]=(*Theta_usAllTreatments.at(tr));
     ListToReturn["Acc_Theta"]=MatAccTheta_k;
     ListToReturn["PostPredictiveY"]=    (*PostPredictive_posteriorAllTreatments.at(tr));
     
     
     std::string treatmentname = "Treatment" + std::to_string(tr + 1);
     
     LISTallTREATMENTS[treatmentname]=ListToReturn;
     
     //
     //InitialVal[treatmentname]=(*STATEMVBPRAllTreatVectOfPtr.at(tr)).GetInitialvalues();
     
   }
   
   
   
   
   if(NprognVar!=0){
     
     LISTallTREATMENTS["Beta"]=*Betas_posterior;
     LISTallTREATMENTS["Betaacc"]=(*BETACLASS).getAcc()/(MCMCTotIter-Burnin);
     
   }
   
   
   //info
   LISTallTREATMENTS["Info"]=Rcpp::List::create(
     Rcpp::_["MCMCinfo"]=Rcpp::NumericVector::create(
                                 Rcpp::_["Burnin"] = Burnin,
                                 Rcpp::_["MCMCFinalSamplesize"] = MCMCFinalSamplesize,
                                 Rcpp::_["Thinning"] =Thinning,
                                 Rcpp::_["Time"] =elapsed ),
     Rcpp::_["InitialValues"]=InitialVal,
     Rcpp::_["VariablesNames"]=ListVariablesNames
   );
   
  // LISTallTREATMENTS["Info:"]= Rcpp::NumericVector::create(
  //   Rcpp::_["Burnin"] = Burnin,
  //   Rcpp::_["MCMCFinalSamplesize"] = MCMCFinalSamplesize,
  //   Rcpp::_["Thinning"] =Thinning,
  //   Rcpp::_["Time"] =elapsed
   
  // );
   
  // LISTallTREATMENTS["initvalues"]=InitialVal;
   
  // LISTallTREATMENTS["Varnames"]=ListVariablesNames;
   
   
   return LISTallTREATMENTS;
   
 };


void test(){};
