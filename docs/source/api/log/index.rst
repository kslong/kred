log
===

.. py:module:: log

.. autoapi-nested-parse::

   This is a general purpose logging utility.  It is intended to record the history
   of associated with processing of the data, mainly but not exclusively of
   errors that are encountered.



Attributes
----------

.. autoapisummary::

   log.LOG_DIR
   log.log_file


Functions
---------

.. autoapisummary::

   log.close_log
   log.get_current_commit
   log.log_message
   log.open_log


Module Contents
---------------

.. py:data:: LOG_DIR
   :value: './Log/'


.. py:function:: close_log()

.. py:function:: get_current_commit()

.. py:data:: log_file
   :value: None


.. py:function:: log_message(message)

.. py:function:: open_log(file_name, reinitialize=False)

