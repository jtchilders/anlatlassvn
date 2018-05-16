#!/usr/bin/env python
import os,sys,optparse,logging
from OpenSSL._util import lib as sslib
logger = logging.getLogger(__name__)

def main():
   ''' simple starter program that can be copied for use when starting a new script. '''
   logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(name)s:%(message)s')

   parser = optparse.OptionParser(description='')
   options,args = parser.parse_args()

   
   manditory_args = [
                     
                  ]

   for man in manditory_args:
      if options.__dict__[man] is None:
         logger.error('Must specify option: ' + man)
         parser.print_help()
         sys.exit(-1)
   

   # get proxy file
   fileName = os.environ['X509_USER_PROXY']
   free_fileName = 0
   if fileName is None:
      fileName = '/tmp/x509up_u%d' % os.getuid()
      free_fileName = 1

   # Find the trusted CA directory 
   CA_dir = os.environ['X509_CERT_DIR']
   if CA_dir is None: CA_dir = "/etc/grid-security/certificates/";

   logger.debug("Testing CA directory %s", CA_dir )
   if not os.path.exists(CA_dir):
      CA_dir = '/.globus/certificates/'
      if not os.path.exists(CA_dir):
         raise Exception(' Trusted certificate directory not found!')

   # Check the file permissions on the proxy
   logger.info('Checking file permissions for ' + fileName)
   statinfo = os.stat(CA_dir)
   if statinfo.st_mode & (statinfo.S_IRWXG | statinfo.S_IRWXO ):
      raise Exception('should be 0600, are currently %04o.'%(statinfo.st_mode & 0xFFF))
   
   # read the file and build the certificate stack, plus the proxy private key 
   certStack = sslib.sk_X509_new_null();
   if certstack is None: return sslib.ERR_get_error()
   result = grid_readProxy( fileName, certStack, pkey )
   if ( result == X509_V_OK ):
      # first, verify the certificate chain 
      result = grid_verifyCert( CA_dir, certStack )
      if ( result == X509_V_OK ):
         # still OK? then match the proxy public and private keys
         result = grid_verifyPrivateKey( certStack, pkey )
         if ( result == X509_V_OK )
            # aaah, nirwana */
            printf( "OK\n" )
         else
            raise Exception( "Verifying private key", "%s\n",ERR_reason_error_string( result ) );
      else
         raise Exception( "Verifying certificate chain", "%s\n",X509_verify_cert_error_string( result ) );
   else
     raise Exception( "Reading proxy", "%s\n", ERR_reason_error_string( result ) );

def grid_readProxy(filename,certstack,pkey):

   sk      = None # STACK_OF(X509_INFO)
   xi      = None # X509_INFO
   err     = 0

   logger.debug("--- Welcome to the grid_readProxy function ---")


   logger.info("Reading file %s", filename )

   certbio = sslib.BIO_new_file( filename, "r" )
   if certbio is None:
      return sslib.ERR_get_error()

   logger.debug("Reading X509_INFO records" )
   sk=sslib.PEM_read_bio_X509(certbio, None, None, None)
   if not sk:
      err = sslib.ERR_get_error()
      logger.error("No X509 records found" )
      sslib.BIO_free(certbio)
      sslib.sk_X509_INFO_free(sk)
      sslib.sk_X509_free(certstack)
      certstack = None
      return err

   logger.debug("Resetting BIO")
   err = sslib.BIO_reset( certbio )
   if ( err != sslib.X509_V_OK ) return err

   logger.debug("Reading Private key" )
   pkey = PEM_read_bio_PrivateKey( certbio, None, grid_X509_empty_callback, None )

   if pkey is None: logger.warning("No private key found." )

   while (sslib.sk_X509_INFO_num(sk)):
      xi=sk_X509_INFO_shift(sk)
      if (xi->x509 != NULL):
         sk_X509_push(*certstack, xi->x509)
         xi->x509=NULL
      X509_INFO_free(xi)
    
   if (!sk_X509_num(*certstack)):
      err = ERR_get_error()
      Error( oper, "No certificates found" )
      BIO_free(certbio)
      sk_X509_INFO_free(sk)
      sk_X509_free(*certstack)
      *certstack = NULL
      return err

   BIO_free(certbio)
   sk_X509_INFO_free(sk)

   return X509_V_OK

def grid_X509_empty_callback(buf,buf_size,verify,cb_tmp):
   if ( buf_size > 0 ) buf = '\0'
   return 0

if __name__ == "__main__":
   main()
