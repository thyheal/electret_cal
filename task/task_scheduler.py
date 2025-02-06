from typing import List, Optional, Dict
from concurrent.futures import ThreadPoolExecutor
import logging
import time
from .task_manager import TaskManager, TaskStatus, Task
from utils.cluster_monitor import ClusterMonitor

class TaskScheduler:
    def __init__(self, max_workers: int = 4, debug: bool = False):
        """初始化任务调度器

        Args:
            max_workers: 最大并行任务数
            debug: 是否为调试模式，调试模式下限制队列任务数为5
        """
        self.max_workers = max_workers
        self.debug = debug
        self.task_manager = TaskManager()
        self.executor = ThreadPoolExecutor(max_workers=max_workers)
        self.running_futures: Dict[str, 'Future'] = {}
        self.cluster_monitor = ClusterMonitor()
        self.logger = logging.getLogger(__name__)

    def schedule_task(self, molecule_name: str, calculation_params: Dict) -> str:
        """调度新任务

        Args:
            molecule_name: 分子名称
            calculation_params: 计算参数

        Returns:
            str: 任务ID
        """
        task_id = self.task_manager.create_task(molecule_name, calculation_params)
        self._try_execute_pending_tasks()
        return task_id

    def _try_execute_pending_tasks(self) -> None:
        """尝试执行待处理的任务"""
        running_tasks = self.task_manager.get_running_tasks()
        if len(running_tasks) >= self.max_workers:
            return

        # 在调试模式下检查集群队列状态
        if self.debug:
            max_cluster_tasks = 5
            while self.cluster_monitor.is_queue_full(max_cluster_tasks):
                self.logger.info("等待集群队列空闲...")
                time.sleep(60)

        available_slots = self.max_workers - len(running_tasks)
        pending_tasks = self.task_manager.get_pending_tasks()

        for task in pending_tasks[:available_slots]:
            self._execute_task(task)

    def _execute_task(self, task: Task) -> None:
        """执行单个任务

        Args:
            task: 要执行的任务
        """
        self.task_manager.update_task_status(task.id, TaskStatus.RUNNING)
        future = self.executor.submit(self._run_calculation, task)
        self.running_futures[task.id] = future
        future.add_done_callback(lambda f: self._handle_task_completion(task.id, f))

    def _run_calculation(self, task: Task) -> None:
        """运行高斯计算

        Args:
            task: 要执行的任务
        """
        try:
            # TODO: 实现与GaussianCalculator的集成
            # 这里将添加实际的计算逻辑
            pass
        except Exception as e:
            self.logger.error(f"Task {task.id} failed: {str(e)}")
            raise

    def _handle_task_completion(self, task_id: str, future: 'Future') -> None:
        """处理任务完成

        Args:
            task_id: 任务ID
            future: 任务的Future对象
        """
        try:
            future.result()  # 如果有异常会在这里抛出
            self.task_manager.update_task_status(task_id, TaskStatus.COMPLETED)
        except Exception as e:
            self.task_manager.update_task_status(
                task_id,
                TaskStatus.FAILED,
                str(e)
            )
        finally:
            if task_id in self.running_futures:
                del self.running_futures[task_id]
            self._try_execute_pending_tasks()

    def get_task_status(self, task_id: str) -> Optional[TaskStatus]:
        """获取任务状态

        Args:
            task_id: 任务ID

        Returns:
            Optional[TaskStatus]: 任务状态
        """
        return self.task_manager.get_task_status(task_id)

    def shutdown(self) -> None:
        """关闭调度器"""
        self.executor.shutdown(wait=True)