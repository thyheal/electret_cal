import subprocess
from typing import Optional

class ClusterMonitor:
    def __init__(self):
        """初始化集群监控器"""
        pass

    def get_running_tasks_count(self, queue_name: str = 'i8cpu') -> int:
        """获取指定队列中正在运行的任务数量

        Args:
            queue_name: 队列名称，默认为'i8cpu'

        Returns:
            int: 正在运行的任务数量
        """
        try:
            output = subprocess.check_output("squeue", shell=True, text=True)
            return output.count(queue_name)
        except subprocess.CalledProcessError:
            return 0

    def is_queue_full(self, max_tasks: int, queue_name: str = 'i8cpu') -> bool:
        """检查指定队列是否已满

        Args:
            max_tasks: 最大允许的任务数量
            queue_name: 队列名称，默认为'i8cpu'

        Returns:
            bool: 如果队列中的任务数量达到或超过最大值则返回True
        """
        return self.get_running_tasks_count(queue_name) >= max_tasks