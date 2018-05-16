#!/usr/bin/env bash
ROOT=/users/hpcusers/argobalsam/gridsecurity
export X509_USER_CERT=$ROOT/$USER/xrootdsrv-cert.pem
export X509_USER_KEY=$ROOT/$USER/xrootdsrv-key.pem
export X509_USER_PROXY=$ROOT/$USER/proxy
export X509_CERT_DIR=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/etc/grid-security-emi/certificates
#export X509_CERT_DIR=$ROOT/certificates
#export X509_CERT_DIR=/etc/grid-security/certificates/
export X509_CA_CERTS=$ROOT/$USER/cacerts.pem
