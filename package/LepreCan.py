#!/usr/bin/env python
# -*- coding: utf-8 -*-
#==========================================================================
# 
#--------------------------------------------------------------------------
# Copyright (c) 2022 Light Conversion, UAB
# All rights reserved.
# www.lightcon.com
#==========================================================================
import numpy
from enum import Enum, IntEnum
import time
import struct

from collections import deque
from threading import Lock        

class LepreCanDevice:
    baseId = 0x100
    
    RegisterResponseTimeout = 1.5 # in [seconds]
    RegisterResponseRetries = 3
    
    mutex = Lock()
    
    def __init__ (self, km, baseId):
        self.km = km
        self.baseId = baseId
    
    def GetRegister(self, register, index, flags=0x00):
        can_id = self.GenerateMessageId(self.baseId, FrameType.GetRegisterCommandFrame)
        
        data = self.GenerateDataFrame(frameType = FrameType.GetRegisterCommandFrame, 
                                  registerAddress = register,
                                  index = index or 0x00, 
                                  flags = 0x00, 
                                  data4bytes=0x00000000)
        
        self.mutex.acquire()
        self.km.send(can_id, data)
        recv = self.km.receive()
        self.mutex.release()
        
        return (recv[3], recv[4:])
    
    def SetRegister(self, register, index, data4bytes, flags=0x00):
        can_id = self.GenerateMessageId(self.baseId, FrameType.SetRegisterCommandFrame)
        
        data = self.GenerateDataFrame(frameType = FrameType.SetRegisterCommandFrame, 
                                  registerAddress = register,
                                  index = index or 0x00, 
                                  flags = 0x00, 
                                  data4bytes=data4bytes)
        
#        usb2can.PrintCanalMsg("send", send_msg)                    
        self.mutex.acquire()
        self.km.send(can_id, data)
        recv = self.km.receive()
        self.mutex.release()
                
        return (recv[3], recv[4:])
        
    def SetRegisterAsync(self, register, index, data4bytes, flags=0x00):
        can_id = self.GenerateMessageId(self.baseId, FrameType.SetRegisterCommandFrame)
        
        data = self.GenerateDataFrame(frameType = FrameType.SetRegisterCommandFrame, 
                                  registerAddress = register,
                                  index = index or 0x00, 
                                  flags = 0x00, 
                                  data4bytes=data4bytes)
        
#        usb2can.PrintCanalMsg("send", send_msg)                    
        self.mutex.acquire()
        self.km.send(can_id, data)
        self.mutex.release()
                
        return (0, 0)

    def float_to_hex(self, f):
        return int(struct.unpack('<I', struct.pack('<f', f))[0])
    
    def hex_to_float(self, f):
        return float(struct.unpack('<f', struct.pack('<I', f))[0])
    
    def IterateCrc8Byte(self, seed, newByte):
        data = numpy.uint8(seed ^ newByte);
        for i in range(8):
            if ((data & 0x80) != 0):
                data <<= 1;
                data ^= 0x07;
            else:
                data <<= 1;    
        return numpy.uint8(data);
    
    def PrepareFrame(self, data):
        indexes = [0, 1, 2, 4, 5, 6, 7]
        
        crc = self.IterateCrc8Byte (0xff, data[0]);    
        
        for index in indexes[1:]:
            crc = self.IterateCrc8Byte (crc, data[index])
            
        data[3] = crc;
        return data
        
    def GetDataString (self, data):
        return ' '.join(["{0:#0{1}x}".format(cell,4)[2:] for cell in data])
    
    
    def GenerateMessageId (self, baseId, frameType):
        return baseId + frameType.value
    
    def GenerateDataFrame (self, frameType, registerAddress, index, flags=0x00, data4bytes=0x00000000):
        data = [0] * 8
                   
        
        if (frameType == FrameType.GetRegisterCommandFrame):
            data = self.PrepareFrame([
                    registerAddress & 0x00ff,
                    (registerAddress & 0xff00) >> 8,
                    index,
                    0x00, # will be replaced by crc
                    data4bytes & 0x000000ff,
                    (data4bytes & 0x0000ff00) >> 8,
                    (data4bytes & 0x00ff0000) >> 16,
                    (data4bytes & 0xff000000) >> 24
                    ]);     
        
        if (frameType == FrameType.SetRegisterCommandFrame):
            data = self.PrepareFrame([
                    registerAddress & 0x00ff,
                    (registerAddress & 0xff00) >> 8,
                    index,
                    0x00, # will be replaced by crc
                    data4bytes & 0x000000ff,
                    (data4bytes & 0x0000ff00) >> 8,
                    (data4bytes & 0x00ff0000) >> 16,
                    (data4bytes & 0xff000000) >> 24
                    ]);   
        return data

def BytesArrayToInt (data):
    result = 0
 
    for i, val in enumerate(data):
        result = result + (val << (i * 8))
                                 
    return result     
        
        
class FrameType(IntEnum):
    BroadcastFrame = 0
    SetRegisterResponseFrame = 1
    GetRegisterResponseFrame = 2
    OutgoingIfuFrame = 3
    SetRegisterCommandFrame = 5
    GetRegisterCommandFrame = 6
    IncomingRawFrame = 7
    Invalid = 255
    
class ResponseStatus(IntEnum) :
    Success =                   0x00,
    UnsupportedCommand =        0x01,
    InvalidRegisterAddress =    0x02,
    InvalidIndex =              0x03,
    TypeError =                 0x04,
    DecryptionFailed =          0x05,
    InsufficientAccessLevel =   0x08,
    UnknownError2 =              0xef,
    UnknownError =              0xff