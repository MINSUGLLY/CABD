#include "simulator.hpp"

namespace CABD{

void Simulator::advance(){
    for(const auto& object : objsList){
        objectAdvance(object);
    }
    for(const auto& multibody : multibodiesList){
        multibodyAdvance(multibody);
    }
}

void Simulator::objectAdvance(std::shared_ptr<Object> object){
    vec3 p;
    mat3 R;
    utils::q2Rp(object->q, R, p);

    vec3 pdot;
    mat3 Rdot;
    utils::q2Rp(object->qdot, Rdot, pdot);

    /* corotated ABD */
    vec3 forceLocal = R.transpose() * (object->forceE + object->massTotal * gravity);
    vec3 torqueLocal = R.transpose() * object->torqueE;
    vec3 velocityLocal = R.transpose() * pdot;
    mat3 rotationLocal = R.transpose() * Rdot;
    vec12 qdotLocal;
    utils::Rp2q(qdotLocal, rotationLocal, velocityLocal);

    vec12 linearGradient = vec12::Zero();

    matJ J;
    utils::getJ(vec3(0.0, 0.0, 0.0), J);
    linearGradient += -J.transpose() * forceLocal * dt * dt - object->Mtet * qdotLocal * dt;
    utils::getJ(vec3(1.0, 0.0, 0.0), J);
    linearGradient += -J.transpose() * vec3(0.0, torqueLocal[2], 0.0) * dt * dt;
    utils::getJ(vec3(0.0, 1.0, 0.0), J);
    linearGradient += -J.transpose() * vec3(0.0, 0.0, torqueLocal[0]) * dt * dt;
    utils::getJ(vec3(0.0, 0.0, 1.0), J);
    linearGradient += -J.transpose() * vec3(torqueLocal[1], 0.0, 0.0) * dt * dt;
    utils::getJ(vec3(0.0, 0.0, 0.0), J);
    linearGradient += J.transpose() * vec3(torqueLocal[1], torqueLocal[2], torqueLocal[0]) * dt * dt;
    
    vec12 result = -object->linearHessianInv * linearGradient;
        
    vec3 dpLocal = result.block<3, 1>(0, 0);
    mat3 dRLocal;
    dRLocal.col(0) = result.block<3, 1>(3, 0) + vec3(1.0, 0.0, 0.0);
    dRLocal.col(1) = result.block<3, 1>(6, 0) + vec3(0.0, 1.0, 0.0);
    dRLocal.col(2) = result.block<3, 1>(9, 0) + vec3(0.0, 0.0, 1.0);

    vec3 dpGlobal = R * dpLocal;

    mat3 polarDecompositionR = utils::mat3PolarDecomposition(R * dRLocal);

    utils::Rp2q(object->qtemp, polarDecompositionR, p + dpGlobal);
    /* END */

    object->qdot = (object->qtemp - object->q) / dt;
    {
        std::lock_guard<std::mutex> lock(mtx);
        object->q = object->qtemp;
    }
}

void Simulator::multibodyAdvance(std::shared_ptr<Multibody> multibody){
    // gradient
    vecX gradient;
    gradient.resize(12 * multibody->linksList.size() + 3 * multibody->jointsList.size());
    gradient.setZero();

    matJ J;
    for(auto& link : multibody->linksList){
        if(!link->fixed){
            vec3 p;
            mat3 R;
            utils::q2Rp(link->q, R, p);

            vec3 pdot;
            mat3 Rdot;
            utils::q2Rp(link->qdot, Rdot, pdot);

            vec3 forceLocal = R.transpose() * (link->forceE + link->massTotal * gravity);
            vec3 torqueLocal = R.transpose() * link->torqueE;
            vec3 velocityLocal = R.transpose() * pdot;
            mat3 rotationLocal = R.transpose() * Rdot;
            vec12 qdotLocal;
            utils::Rp2q(qdotLocal, rotationLocal, velocityLocal);

            utils::getJ(vec3(0.0, 0.0, 0.0), J);
            gradient.block<12, 1>(12 * link->idx, 0) += -J.transpose() * forceLocal * dt * dt - link->Mtet * qdotLocal * dt;
        }
        else{
            vec3 p;
            mat3 R;
            utils::q2Rp(link->q, R, p);

            link->pt += link->dp;

            vec12 qtlocal;
            utils::Rp2q(qtlocal, R.transpose() * link->Rt, R.transpose() * (link->pt - p));

            vec12 qo;
            utils::Rp2q(qo, mat3::Identity(), vec3::Zero());

            gradient.block<12, 1>(12 * link->idx, 0) += -link->Mtet * (qtlocal - qo);
        }
    }

    if(multibody->jointsList.size() != 0){
        // set constraint
        matX matC;
        matC.resize(3 * multibody->jointsList.size(), 12 * multibody->linksList.size());
        matC.setZero();
        int rowIdx = 0;
        for(const auto& joint : multibody->jointsList){
            vec3 p0;
            mat3 R0;
            utils::q2Rp(multibody->linksList[joint->idx0]->q, R0, p0);

            vec3 p1;
            mat3 R1;
            utils::q2Rp(multibody->linksList[joint->idx1]->q, R1, p1);

            matJ JR0;
            matJ JR1;
            matJ J0;
            matJ J1;

            utils::getJR(joint->pos0, R0, JR0);
            matC.block<3, 12>(rowIdx, 12 * multibody->linksList[joint->idx0]->idx) = JR0;
            utils::getJR(joint->pos1, R1, JR1);
            matC.block<3, 12>(rowIdx, 12 * multibody->linksList[joint->idx1]->idx) = -JR1;
            
            utils::getJ(joint->pos0, J0);
            utils::getJ(joint->pos1, J1);
            gradient.block<3, 1>(12 * multibody->linksList.size() + rowIdx, 0) = J0 * multibody->linksList[joint->idx0]->q - J1 * multibody->linksList[joint->idx1]->q;
        
            rowIdx += 3;
        }

        matX matS;
        matS.resize(3 * multibody->jointsList.size(), 3 * multibody->jointsList.size());
        matS.setZero();
        vecX CHinvf;
        CHinvf.resize(3 * multibody->jointsList.size());
        CHinvf.setZero();
        for(const auto& link : multibody->linksList){
            matX Cn = matC.block(0, 12 * link->idx, 3 * multibody->jointsList.size(), 12);
            
            for(int i = 0; i < link->begList.size(); i++){
                for(int j = i; j < link->begList.size(); j++){
                    matX temp = Cn.block(link->begList[i], 0, 3, 12) * link->linearHessianInv * Cn.block(link->begList[j], 0, 3, 12).transpose();
                    matS.block(link->begList[i], link->begList[j], 3, 3) += temp;
                    if(i != j){
                        matS.block(link->begList[j], link->begList[i], 3, 3) += temp.transpose();
                    }
                }
                CHinvf.block(link->begList[i], 0, 3, 1) += Cn.block(link->begList[i], 0, 3, 12) * link->linearHessianInv * gradient.block<12, 1>(12 * link->idx, 0);
            }
        }
        vecX g_CHinvf = gradient.block(12 * multibody->linksList.size(), 0, 3 * multibody->jointsList.size(), 1) - CHinvf;

        /* LLT */
        Eigen::LLT<Eigen::MatrixXd> llt(matS);
        vecX lam = llt.solve(g_CHinvf);
        /* LLT END*/

        vecX Clam = matC.transpose() * lam;

        for(auto& link : multibody->linksList){
            vec3 p;
            mat3 R;
            utils::q2Rp(link->q, R, p);
            vec12 resulttemp = -link->linearHessianInv * (gradient.block<12, 1>(12 * link->idx, 0) + Clam.block<12, 1>(12 * link->idx, 0));

            vec3 dpLocal = resulttemp.block<3, 1>(0, 0);
            mat3 dRLocal;
            dRLocal.col(0) = resulttemp.block<3, 1>(3, 0) + vec3(1.0, 0.0, 0.0);
            dRLocal.col(1) = resulttemp.block<3, 1>(6, 0) + vec3(0.0, 1.0, 0.0);
            dRLocal.col(2) = resulttemp.block<3, 1>(9, 0) + vec3(0.0, 0.0, 1.0);

            vec3 dpGlobal = R * dpLocal;

            mat3 polarDecompositionR = utils::mat3PolarDecomposition(R * dRLocal); 
            
            utils::Rp2q(link->qtemp, polarDecompositionR, p + dpGlobal);
            link->qdot = (link->qtemp - link->q) / dt;
            {
                std::lock_guard<std::mutex> lock(mtx);
                link->q = link->qtemp;
            }
        }
    }
    else{
        for(auto& link : multibody->linksList){
            vec3 p;
            mat3 R;
            utils::q2Rp(link->q, R, p);
            vec12 resulttemp = -link->linearHessianInv * gradient.block<12, 1>(12 * link->idx, 0);

            vec3 dpLocal = resulttemp.block<3, 1>(0, 0);
            mat3 dRLocal;
            dRLocal.col(0) = resulttemp.block<3, 1>(3, 0) + vec3(1.0, 0.0, 0.0);
            dRLocal.col(1) = resulttemp.block<3, 1>(6, 0) + vec3(0.0, 1.0, 0.0);
            dRLocal.col(2) = resulttemp.block<3, 1>(9, 0) + vec3(0.0, 0.0, 1.0);

            vec3 dpGlobal = R * dpLocal;

            mat3 polarDecompositionR = utils::mat3PolarDecomposition(R * dRLocal); 
            
            utils::Rp2q(link->qtemp, polarDecompositionR, p + dpGlobal);
            link->qdot = (link->qtemp - link->q) / dt;
            {
                std::lock_guard<std::mutex> lock(mtx);
                link->q = link->qtemp;
            }
        }
    }
}



}