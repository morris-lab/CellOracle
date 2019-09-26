# -*- coding: utf-8 -*-



import io
import logging
import os


from datetime import datetime


class makelog():
    """
    This is a class for making log.

    """
    def __init__(self, file_name=None, directory=None):

        name = datetime.now().ctime()

        # ログの出力名を設定
        logger = logging.Logger('LoggingTest1')

        # ログレベルの設定
        logger.setLevel(10)

        # ログのファイル出力先を設定

        if not file_name is None:
            name = file_name + " " + name

        file_path = f'log {name}.log'

        if not directory is None:
            os.makedirs(directory, exist_ok=True)
            file_path = os.path.join(directory, file_path)

        fh = logging.FileHandler(file_path)
        logger.addHandler(fh)

        # ログのコンソール出力の設定
        sh = logging.StreamHandler()
        logger.addHandler(sh)

        # ログの出力形式の設定
        formatter = logging.Formatter('%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
        fh.setFormatter(formatter)
        sh.setFormatter(formatter)

        self.logger = logger

    def info(self, comment):
        """
        Add comment into the log file.

        Args:
            comment (str): comment.
            
        """
        self.logger.log(20, comment)
