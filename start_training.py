import asyncio, websockets, json

async def start_training():
    uri = "ws://localhost:8081/cloud/ws"
    async with websockets.connect(uri) as ws:
        # Запрос на начало actual training (после подтверждения данных)
        msg = {
            "operation": "initialize_cloud_training",  # или "perform_model_evaluation"
            "data": {
                "start_date": "2012-10-12",
                "end_date": "2012-10-13",
                "is_cache_active": True,
                "genetic_evaluation_strategy": "elitism",
                "model_type": "lstm"
            }
        }
        await ws.send(json.dumps(msg))
        print("Ответ:", await ws.recv())

asyncio.run(start_training())
