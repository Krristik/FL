import asyncio
import websockets
import json

async def send_training_request():
    uri = "ws://localhost:8081/cloud/ws"
    async with websockets.connect(uri) as websocket:
        # Сообщение для запуска pretraining
        message = {
            "operation": "initialize_cloud_pretraining",
            "data": {
                "start_date": "2023-01-01",
                "end_date": "2023-01-31",
                "is_cache_active": True,
                "genetic_evaluation_strategy": "elitism",
                "model_type": "lstm"
            }
        }
        await websocket.send(json.dumps(message))
        response = await websocket.recv()
        print(f"Ответ сервера: {response}")

asyncio.run(send_training_request())
