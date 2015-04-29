
''' This file is really just a convenience in that it makes use of two different custom parsers that could have been individually specified '''

import KnownGenesParser, SNP2Signal

class AnnotationParser(KnownGenesParser.AnnotationParser):
  pass

class SignalParser(SNP2Signal.SignalParser):
  pass


