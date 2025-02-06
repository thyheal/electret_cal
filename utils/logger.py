import logging
import os
from datetime import datetime
from typing import Optional

class LoggerConfig:
    def __init__(self,
                 log_level: int = logging.INFO,
                 log_file: Optional[str] = None,
                 console_output: bool = True,
                 log_format: str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'):
        """初始化日志配置

        Args:
            log_level: 日志级别
            log_file: 日志文件路径，如果为None则不输出到文件
            console_output: 是否输出到控制台
            log_format: 日志格式
        """
        self.log_level = log_level
        self.log_file = log_file
        self.console_output = console_output
        self.log_format = log_format

def setup_logger(name: str, config: LoggerConfig) -> logging.Logger:
    """设置日志记录器

    Args:
        name: 日志记录器名称
        config: 日志配置

    Returns:
        logging.Logger: 配置好的日志记录器
    """
    logger = logging.getLogger(name)
    logger.setLevel(config.log_level)

    formatter = logging.Formatter(config.log_format)

    # 清除现有的处理器
    logger.handlers.clear()

    # 添加文件处理器
    if config.log_file:
        log_dir = os.path.dirname(config.log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
        file_handler = logging.FileHandler(config.log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # 添加控制台处理器
    if config.console_output:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    return logger

def get_default_logger(name: str) -> logging.Logger:
    """获取默认配置的日志记录器

    Args:
        name: 日志记录器名称

    Returns:
        logging.Logger: 默认配置的日志记录器
    """
    log_dir = "logs"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    log_file = os.path.join(
        log_dir,
        f"{name}_{datetime.now().strftime('%Y%m%d')}.log"
    )

    config = LoggerConfig(
        log_level=logging.INFO,
        log_file=log_file,
        console_output=True
    )

    return setup_logger(name, config)