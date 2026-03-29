import os
from typing import Optional, List

from shared.fed_node.fed_node import FedNode, ParentFedNode, ChildFedNode, FedNodeType


class NodeState:
    current_node: Optional[FedNode] = None
    _auto_init_attempted: bool = False  # Флаг, чтобы не пытаться инициализировать бесконечно

    @staticmethod
    def initialize_node(node: FedNode) -> None:
        """Ручная инициализация узла (оригинальный метод, не меняем)"""
        if NodeState.current_node is not None:
            raise ValueError("Current working node is already initialized.")
        NodeState.current_node = node

    @staticmethod
    def get_current_node() -> Optional[FedNode]:
        """Возвращает текущий узел, при необходимости авто-инициализируя из env"""
        # Если узел уже задан — возвращаем его
        if NodeState.current_node is not None:
            return NodeState.current_node

        # Если ещё не пробовали авто-инициализацию — пробуем
        if not NodeState._auto_init_attempted:
            NodeState._auto_init_attempted = True
            NodeState._try_auto_initialize_from_env()

        return NodeState.current_node

    @staticmethod
    def _try_auto_initialize_from_env() -> None:
        """Пытается создать FedNode на основе переменных окружения (для симуляции)"""
        node_type = os.getenv("NODE_TYPE")  # "cloud", "fog", "edge"
        node_id = os.getenv("NODE_ID", "default")  # "edge_1", "edge_2", etc.

        if not node_type:
            # Нет переменной — не можем авто-инициализировать, возвращаем None
            return

        # Порт всегда 8081 (настроен в Dockerfile)
        port = 8081

        # Создаём базовую конфигурацию для симуляции
        # Используем имена сервисов Docker вместо IP-адресов
        if node_type == "cloud":
            # Cloud: родитель — None, дети — fog-узлы
            node = ParentFedNode(
                node_id="cloud_node_001",
                name="cloud_app",
                fed_node_type=FedNodeType.CLOUD_NODE,
                ip_address="rabbitmq",  # RabbitMQ доступен по имени сервиса
                port=port
            )
            # Добавляем Fog как дочерний узел
            fog_child = ChildFedNode(
                node_id="fog_node_001",
                name="fog_app",
                fed_node_type=FedNodeType.FOG_NODE,
                ip_address="fog_app",
                port=port
            )
            node.add_child_node(fog_child)

        elif node_type == "fog":
            # Fog: родитель — cloud, дети — edge-узлы
            node = ParentFedNode(
                node_id="fog_node_001",
                name="fog_app",
                fed_node_type=FedNodeType.FOG_NODE,
                ip_address="rabbitmq",
                port=port
            )
            # Добавляем Cloud как родителя
            cloud_parent = ParentFedNode(
                node_id="cloud_node_001",
                name="cloud_app",
                fed_node_type=FedNodeType.CLOUD_NODE,
                ip_address="cloud_app",
                port=port
            )
            node.set_parent_node(cloud_parent)

            # Добавляем Edge узлы как детей
            edge_children = [
                ChildFedNode(
                    node_id="edge_node_edge_1",
                    name="edge_app_1",
                    fed_node_type=FedNodeType.EDGE_NODE,
                    ip_address="edge_app_1",
                    port=port
                ),
                ChildFedNode(
                    node_id="edge_node_edge_2",
                    name="edge_app_2",
                    fed_node_type=FedNodeType.EDGE_NODE,
                    ip_address="edge_app_2",
                    port=port
                ),
            ]
            node.add_child_nodes(edge_children)

        elif node_type == "edge":
            # Edge: родитель — fog
            node = ChildFedNode(
                node_id=f"edge_node_{node_id}" if node_id != "default" else "edge_node_001",
                name=f"edge_app_{node_id}" if node_id != "default" else "edge_app",
                fed_node_type=FedNodeType.EDGE_NODE,
                ip_address="rabbitmq",
                port=port
            )
            # Добавляем Fog как родителя
            fog_parent = ParentFedNode(
                node_id="fog_node_001",
                name="fog_app",
                fed_node_type=FedNodeType.FOG_NODE,
                ip_address="fog_app",
                port=port
            )
            node.set_parent_node(fog_parent)
        else:
            # Неизвестный тип узла
            return

        # Инициализируем статическое поле
        NodeState.current_node = node

    @staticmethod
    def reset_node() -> None:
        """Сброс узла (для тестов)"""
        NodeState.current_node = None
        NodeState._auto_init_attempted = False