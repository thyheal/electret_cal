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

    def cancel_all_jobs(self) -> None:
        """取消所有正在排队的任务"""
        try:
            # 获取所有任务的 JOBID
            output = subprocess.check_output("squeue", shell=True, text=True)
            # 跳过表头行，按行分割
            lines = output.strip().split('\n')[1:]
            
            # 提取每个任务的 JOBID 并取消
            for line in lines:
                if line.strip():
                    job_id = line.split()[0]
                    subprocess.run(f"scancel {job_id}", shell=True)
                    print(f"已取消任务: {job_id}")
                    
        except subprocess.CalledProcessError as e:
            print(f"取消任务时出错: {str(e)}")
