from typing import List, Dict, Optional
from dataclasses import dataclass
from datetime import datetime
import logging
from enum import Enum

class TaskStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"

@dataclass
class Task:
    id: str
    molecule_name: str
    calculation_params: Dict
    status: TaskStatus
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None

class TaskManager:
    def __init__(self):
        self.tasks: Dict[str, Task] = {}
        self.logger = logging.getLogger(__name__)

    def create_task(self, molecule_name: str, calculation_params: Dict) -> str:
        """创建新的计算任务

        Args:
            molecule_name: 分子名称
            calculation_params: 计算参数字典

        Returns:
            str: 任务ID
        """
        task_id = f"{molecule_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        task = Task(
            id=task_id,
            molecule_name=molecule_name,
            calculation_params=calculation_params,
            status=TaskStatus.PENDING,
            created_at=datetime.now()
        )
        self.tasks[task_id] = task
        self.logger.info(f"Created task {task_id} for molecule {molecule_name}")
        return task_id

    def get_task(self, task_id: str) -> Optional[Task]:
        """获取任务信息

        Args:
            task_id: 任务ID

        Returns:
            Optional[Task]: 任务对象，如果不存在则返回None
        """
        return self.tasks.get(task_id)

    def update_task_status(self, task_id: str, status: TaskStatus, error_message: Optional[str] = None) -> None:
        """更新任务状态

        Args:
            task_id: 任务ID
            status: 新的任务状态
            error_message: 错误信息（如果有）
        """
        if task_id not in self.tasks:
            raise ValueError(f"Task {task_id} not found")

        task = self.tasks[task_id]
        task.status = status

        if status == TaskStatus.RUNNING and not task.started_at:
            task.started_at = datetime.now()
        elif status in [TaskStatus.COMPLETED, TaskStatus.FAILED]:
            task.completed_at = datetime.now()
            if error_message:
                task.error_message = error_message

        self.logger.info(f"Updated task {task_id} status to {status.value}")

    def get_pending_tasks(self) -> List[Task]:
        """获取所有待处理的任务

        Returns:
            List[Task]: 待处理任务列表
        """
        return [task for task in self.tasks.values() if task.status == TaskStatus.PENDING]

    def get_running_tasks(self) -> List[Task]:
        """获取所有正在运行的任务

        Returns:
            List[Task]: 运行中任务列表
        """
        return [task for task in self.tasks.values() if task.status == TaskStatus.RUNNING]

    def get_task_status(self, task_id: str) -> Optional[TaskStatus]:
        """获取任务状态

        Args:
            task_id: 任务ID

        Returns:
            Optional[TaskStatus]: 任务状态，如果任务不存在则返回None
        """
        task = self.get_task(task_id)
        return task.status if task else None