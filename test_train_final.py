import asyncio, websockets, json

async def send_training_request():
    uri = "ws://localhost:8081/cloud/ws"
    async with websockets.connect(uri) as ws:
        # Используем дату из edge_data_1 (2012-10-12)
        msg = {
            "operation": "initialize_cloud_pretraining",
            "data": {
                "start_date": "2012-10-12",
                "end_date": "2012-10-12",
                "is_cache_active": True,
                "genetic_evaluation_strategy": "elitism",
                "model_type": "lstm"
            }
        }
        await ws.send(json.dumps(msg))
        print("Ответ сервера:", await ws.recv())

asyncio.run(send_training_request())
