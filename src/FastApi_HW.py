from fastapi import FastAPI,status,HTTPException
from models import User

app = FastAPI()


@app.get("/")
def read_root():
    return {"Message": "Welcome to FastAPI"}


users_db = [{"user_id":1,"name":"Maria","location":"London"},
            {"user_id":2,"name":"Alina","location":"Spain"},
            {"user_id":3,"name":"Sasha","location":"London"}]

@app.get("/users")
def retrive_user():
    """Retrieve all users"""
    return users_db

@app.get("/items/{item_id}")
def retrive_item_id(item_id: int):
    return {"item_id": item_id}


@app.get("/users",tags=["users"])
def retrive_users():
    """Retrieve all users"""
    return users_db

@app.post("/add", status_code =status.HTTP_201_CREATED,tags=["users"],summary="create user",response_description="user is created")
def add_user(user:User):
    users_db.append(user)
    return user

@app.put("/add/{user_id}",tags=["users"])
def update_user(user_id:int, updated_user:User):
    for index,user in enumerate(users_db):
        if user['user_id'] == user_id:
            users_db[index] = updated_user.dict()
            return updated_user.dict()
    raise HTTPException(status_code = 404, detail="User is not found")

@app.delete("/users/{user_id}")
def delete_user(user_id:int):
    for index,user in enumerate(users_db):
        if user['user_id'] == user_id:
            deleted_user = users_db.pop(index)
            return deleted_user
    raise HTTPException(status_code = 404, detail="User is not found")
