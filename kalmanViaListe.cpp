// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
const double M_LOG2PI = 1.8378770664093454835606594728112352797227949472756;

// [[Rcpp::export]]
List kalman(
            const List Tl,
            const List Ql,
            const List Q1,
            const List Zl,
            const List Hl,
            const List Rl,
            const arma::sp_mat a1,
            const arma::sp_mat P1,
            const arma::mat Yt,
            const List Smo = "no"
) 
  {
  
  arma::sp_mat T = Rcpp::as<arma::sp_mat>(Tl[0]);
  arma::sp_mat H = Rcpp::as<arma::sp_mat>(Hl[0]);
  arma::sp_mat Q = Rcpp::as<arma::sp_mat>(Ql[0]);
  arma::sp_mat Q1r = Rcpp::as<arma::sp_mat>(Q1[0]);
  arma::sp_mat Z;
  arma::sp_mat R;
  const arma::uword n = Yt.n_cols; // n. of observations
  const arma::uword m = Yt.n_rows; // n. of observed variables
  const arma::uword p = T.n_rows;  // n. of state variables
  
  const arma::uword rQ = Q1r.n_rows;// righe matrice Q
  
  double loglik=0;
  arma::sp_mat att(p, n);
  arma::sp_mat at(p, n+1);
  arma::sp_mat Ythat(m, n);
  arma::sp_mat v(m, n);
  arma::cube Ptt(p, p, n);
  arma::cube Pt(p, p, n+1);
  arma::cube F(m, m, n);
  arma::cube K(p, m, n);
  arma::sp_mat PZt(p, m);
  const arma::sp_mat Tt = T.t();
  arma::sp_mat iF(m, m);
  arma::uvec notna;
  arma::colvec v0;
  arma::sp_mat PZ0t;
  at.col(0) = a1;
  Pt.slice(0) = P1;
  
  // smoother
  arma::sp_mat r(p, n);
  arma::sp_mat alphahat(p,n);
  arma::cube N(p, p, n);
  arma::cube Vt(p, p, n);
  arma::cube Lt(p, p, n);
  arma::sp_mat ZiF(p, p);
  arma::cube Kt(p, m, n);
  
  arma::colvec r0;
  arma::sp_mat N0(p,p);
  
  // smoother disurbance
  arma::sp_mat ut(m,n);
  arma::sp_mat epshat(m,n);
  arma::sp_mat etahat(rQ,n);
  arma::cube D(m,m,n);
  arma::sp_mat V_eps(m,n);
  arma::cube V_eta(rQ,rQ,n+1);
  
    for (arma::uword t=0, k=0; t < n && k < n; ++t, ++k) {
      if (Tl.size() > 1){
        T = Rcpp::as<arma::sp_mat>(Tl[k]);
      }
      else {
        T = Rcpp::as<arma::sp_mat>(Tl[0]);
      }
      if (Ql.size() > 1){
        Q = Rcpp::as<arma::sp_mat>(Ql[k]);
      }
      else {
        Q = Rcpp::as<arma::sp_mat>(Ql[0]);
      }
      if (Zl.size() > 1){
        Z = Rcpp::as<arma::sp_mat>(Zl[k]);
      }
      else {
        Z = Rcpp::as<arma::sp_mat>(Zl[0]);
      }
      if (Hl.size() > 1){
        H = Rcpp::as<arma::sp_mat>(Hl[k]);
      }
      else {
        H = Rcpp::as<arma::sp_mat>(Hl[0]);
      }
      notna = find_finite(Yt.col(t));
      if (notna.n_elem == m) { // no missing obs in y
        // innovation
        Ythat.col(t) = Z*at.col(t);
        v.col(t) = Yt.col(t) - Ythat.col(t);
        PZt = Pt.slice(t)*Z.t();
        F.slice(t) = Z*PZt + H;
        iF = inv_sympd(F.slice(t));
        loglik += (log_det_sympd(F.slice(t)) +
          v.col(t).t()*iF*v.col(t) +
          m*M_LOG2PI).eval().at(0,0);
        // update
        K.slice(t) = PZt * iF;
        att.col(t) = at.col(t) + K.slice(t)*v.col(t);
        Ptt.slice(t) = Pt.slice(t) - PZt*K.slice(t).t();
        // prediction
        at.col(t+1) = T*att.col(t);
        Pt.slice(t+1) = T*Ptt.slice(t)*Tt + Q;
        
        continue;
      }
      if (notna.n_elem == 0) { // all missing obs in y
        // innovation
        Ythat.col(t) = Z*at.col(t);
        v.col(t) = Yt.col(t) - Ythat.col(t);
        PZt = Pt.slice(t)*Z.t();
        F.slice(t) = Z*PZt + H;
        // update
        att.col(t) = at.col(t);
        Ptt.slice(t) = Pt.slice(t);
        // prediction
        at.col(t+1) = T*att.col(t);
        Pt.slice(t+1) = T*Ptt.slice(t)*Tt + Q;
        
        continue;
      }
      // partially missing values in Y
      // innovation
      Ythat.col(t) = Z*at.col(t);
      v.col(t) = Yt.col(t) - Ythat.col(t);
      v0 = arma::colvec(v.col(t))(notna);
      PZt = Pt.slice(t)*Z.t();
      PZ0t = PZt.cols(notna);
      F.slice(t) = Z*PZt + H;
      iF = inv_sympd(F.slice(t)(notna, notna));
      loglik += (log_det_sympd(F.slice(t)(notna, notna)) +
        v0.t()*iF*v0).eval()(0,0) + notna.n_elem*M_LOG2PI;
      // update
      K.slice(t).cols(notna) = PZ0t * iF;
      att.col(t) = at.col(t) + K.slice(t).cols(notna)*v0;
      Ptt.slice(t) = Pt.slice(t) - PZ0t*K.slice(t).cols(notna).t();
      // prediction
      at.col(t+1) = T*att.col(t);
      Pt.slice(t+1) = T*Ptt.slice(t)*Tt + Q;
      
      }
    
    if (Smo.size() == 1){
    if (Rcpp::as<std::string>(Smo[0]) == "Smoother") {
      for(arma::uword t=n-1, p=n-1; t>=1 && p>=1; --t, p--){
        if (Zl.size() > 1){
          Z = Rcpp::as<arma::sp_mat>(Zl[p]);
        }
        else {
          Z = Rcpp::as<arma::sp_mat>(Zl[0]);
        }
        if (Tl.size() > 1){
          T = Rcpp::as<arma::sp_mat>(Tl[p]);
        }
        else {
          T = Rcpp::as<arma::sp_mat>(Tl[0]);
        }
        iF = inv_sympd(F.slice(t));
        Kt.slice(t) = T*Pt.slice(t)*Z.t()*iF;
        Lt.slice(t) = T - Kt.slice(t)*Z;
        ZiF = Z.t()*iF;
        r.col(t-1) = ZiF*v.col(t) + Lt.slice(t).t()*r.col(t); 
        N.slice(t-1) = ZiF*Z + Lt.slice(t).t()*N.slice(t)*Lt.slice(t);
        alphahat.col(t) = at.col(t) + Pt.slice(t)*r.col(t-1);   
        Vt.slice(t) = Pt.slice(t) - Pt.slice(t)*N.slice(t-1)*Pt.slice(t); 
      }
      Z = Rcpp::as<arma::sp_mat>(Zl[0]);
      T = Rcpp::as<arma::sp_mat>(Tl[0]);
      iF = inv_sympd(F.slice(0));
      Kt.slice(0) = T*Pt.slice(0)*Z.t()*iF;
      Lt.slice(0) = T - Kt.slice(0)*Z;
      ZiF = Z.t()*iF;
      r0 = ZiF*v.col(0) + Lt.slice(0).t()*r.col(0); 
      N0 = ZiF*Z + Lt.slice(0).t()*N.slice(0)*Lt.slice(0); 
      alphahat.col(0) = at.col(0) + Pt.slice(0)*r0;    
      Vt.slice(0) = Pt.slice(0) - Pt.slice(0)*N0*Pt.slice(0);
    }

    if (Rcpp::as<std::string>(Smo[0]) == "Disturbance"){
      for(arma::uword t=n-1, p=n-1; t>=1 && p>=1; --t, --p){
        if (Zl.size() > 1){
          Z = Rcpp::as<arma::sp_mat>(Zl[p]);
        }
        else {
          Z = Rcpp::as<arma::sp_mat>(Zl[0]);
        }
        if (Tl.size() > 1){
          T = Rcpp::as<arma::sp_mat>(Tl[p]);
        }
        else {
          T = Rcpp::as<arma::sp_mat>(Tl[0]);
        }
        if (Hl.size() > 1){
          H = Rcpp::as<arma::sp_mat>(Hl[p]);
        }
        else {
          H = Rcpp::as<arma::sp_mat>(Hl[0]);
        }
        if (Q1.size() > 1){
          Q = Rcpp::as<arma::sp_mat>(Q1[p]);
        }
        else {
          Q = Rcpp::as<arma::sp_mat>(Q1[0]);
        }
        if (Rl.size() > 1){
          R = Rcpp::as<arma::sp_mat>(Rl[p]);
        }
        else {
          R = Rcpp::as<arma::sp_mat>(Rl[0]);
        }
      iF = inv_sympd(F.slice(t));
      Kt.slice(t) = T*Pt.slice(t)*Z.t()*iF;
      ut.col(t) = iF*v.col(t) - Kt.slice(t).t()*r.col(t);
      r.col(t-1) = Z.t()*ut.col(t) + T.t()*r.col(t);
      epshat.col(t) = H*ut.col(t);
      etahat.col(t) = Q*R.t()*r.col(t);
      D.slice(t) = iF + Kt.slice(t).t()*N.slice(t)*Kt.slice(t);
      N.slice(t-1) =  Z.t()*D.slice(t)*Z + T.t()*N.slice(t)*T - Z.t()*Kt.slice(t).t()*N.slice(t)*T - T.t()*N.slice(t)*Kt.slice(t)*Z;
      V_eps.col(t) = (H - H*D.slice(t)*H).diag();
      V_eta.slice(t) = Q - Q*R.t()*N.slice(t)*R*Q;
    }
      Z = Rcpp::as<arma::sp_mat>(Zl[0]);
      T = Rcpp::as<arma::sp_mat>(Tl[0]);
      H = Rcpp::as<arma::sp_mat>(Hl[0]);
      Q = Rcpp::as<arma::sp_mat>(Q1[0]);
      R = Rcpp::as<arma::sp_mat>(Rl[0]);
      iF = inv_sympd(F.slice(0));
      Kt.slice(0) = T*Pt.slice(0)*Z.t()*iF;
      ut.col(0) = iF*v.col(0) - Kt.slice(0).t()*r.col(0);
      epshat.col(0) = H*ut.col(0);
      etahat.col(0) = Q*R.t()*r.col(0);
      D.slice(0) = iF + Kt.slice(0).t()*N.slice(0)*Kt.slice(0);
      V_eps.col(0) = (H - H*D.slice(0)*H).diag();
      V_eta.slice(0) = Q - Q*R.t()*N.slice(0)*R*Q;
    }    
    }
    if (Smo.size() == 2){
      if (((Rcpp::as<std::string>(Smo[0]) == "Smoother") and (Rcpp::as<std::string>(Smo[1]) == "Disturbance"))
            or ((Rcpp::as<std::string>(Smo[0]) == "Disturbance") and (Rcpp::as<std::string>(Smo[1]) == "Smoother"))){
        for(arma::uword t=n-1, p=n-1; t>=1 && p>=1; --t, p--){
          if (Zl.size() > 1){
            Z = Rcpp::as<arma::sp_mat>(Zl[p]);
          }
          else {
            Z = Rcpp::as<arma::sp_mat>(Zl[0]);
          }
          if (Tl.size() > 1){
            T = Rcpp::as<arma::sp_mat>(Tl[p]);
          }
          else {
            T = Rcpp::as<arma::sp_mat>(Tl[0]);
          }
          if (Hl.size() > 1){
            H = Rcpp::as<arma::sp_mat>(Hl[p]);
          }
          else {
            H = Rcpp::as<arma::sp_mat>(Hl[0]);
          }
          if (Q1.size() > 1){
            Q = Rcpp::as<arma::sp_mat>(Q1[p]);
          }
          else {
            Q = Rcpp::as<arma::sp_mat>(Q1[0]);
          }
          if (Rl.size() > 1){
            R = Rcpp::as<arma::sp_mat>(Rl[p]);
          }
          else {
            R = Rcpp::as<arma::sp_mat>(Rl[0]);
          }
          iF = inv_sympd(F.slice(t));
          Kt.slice(t) = T*Pt.slice(t)*Z.t()*iF;
          Lt.slice(t) = T - Kt.slice(t)*Z;
          ZiF = Z.t()*iF;
          r.col(t-1) = ZiF*v.col(t) + Lt.slice(t).t()*r.col(t); 
          N.slice(t-1) = ZiF*Z + Lt.slice(t).t()*N.slice(t)*Lt.slice(t);
          
          alphahat.col(t) = at.col(t) + Pt.slice(t)*r.col(t-1);   
          Vt.slice(t) = Pt.slice(t) - Pt.slice(t)*N.slice(t-1)*Pt.slice(t);
          
          ut.col(t) = iF*v.col(t) - Kt.slice(t).t()*r.col(t);
          epshat.col(t) = H*ut.col(t);
          etahat.col(t) = Q*R.t()*r.col(t);
          D.slice(t) = iF + Kt.slice(t).t()*N.slice(t)*Kt.slice(t);
          V_eps.col(t) = (H - H*D.slice(t)*H).diag();
          V_eta.slice(t) = Q - Q*R.t()*N.slice(t)*R*Q;
        }
        Z = Rcpp::as<arma::sp_mat>(Zl[0]);
        T = Rcpp::as<arma::sp_mat>(Tl[0]);
        H = Rcpp::as<arma::sp_mat>(Hl[0]);
        Q = Rcpp::as<arma::sp_mat>(Q1[0]);
        R = Rcpp::as<arma::sp_mat>(Rl[0]);
        
        iF = inv_sympd(F.slice(0));
        Kt.slice(0) = T*Pt.slice(0)*Z.t()*iF;
        Lt.slice(0) = T - Kt.slice(0)*Z;
        ZiF = Z.t()*iF;
        r0 = ZiF*v.col(0) + Lt.slice(0).t()*r.col(0); 
        N0 = ZiF*Z + Lt.slice(0).t()*N.slice(0)*Lt.slice(0);
        
        alphahat.col(0) = at.col(0) + Pt.slice(0)*r0;    
        Vt.slice(0) = Pt.slice(0) - Pt.slice(0)*N0*Pt.slice(0);
        
        ut.col(0) = iF*v.col(0) - Kt.slice(0).t()*r.col(0);
        epshat.col(0) = H*ut.col(0);
        etahat.col(0) = Q*R.t()*r.col(0);
        D.slice(0) = iF + Kt.slice(0).t()*N.slice(0)*Kt.slice(0);
        
        V_eps.col(0) = (H - H*D.slice(0)*H).diag();
        V_eta.slice(0) = Q - Q*R.t()*N.slice(0)*R*Q;
      }
    }
    
    return List::create(Named("at") = at,
                      Named("Pt") = Pt,
                      Named("att") = att,
                      Named("Ptt") = Ptt,
                      Named("Yhat") = Ythat,
                      Named("v") = v,
                      Named("F") = F,
                      Named("K") = K,
                      Named("loglik") = -0.5*loglik,
                      Named("Vt") = Vt,
                      Named("V_eta") = V_eta,
                      Named("alphahat") = alphahat,
                      Named("V_eps") = V_eps,
                      Named("epshat") = epshat,
                      Named("etahat") = etahat
  );
}