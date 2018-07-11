# coding: utf-8

import logging

logger = logging.getLogger('root')

class ConfigFileError(Exception):
  def __init__(self,*args,**kwargs):
    logger.error(str(*args))
    super().__init__(*args,**kwargs)

class FilteringError(Exception):
  def __init__(self,*args,**kwargs):
    logger.error(str(*args))
    super().__init__(*args,**kwargs)

class STANError(Exception):
  def __init__(self,*args,**kwargs):
    logger.error(str(*args))
    super().__init__(*args,**kwargs)